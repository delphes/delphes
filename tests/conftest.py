import contextlib
import math
import subprocess
import sys
from pathlib import Path

import pytest
import ROOT

from delphes.dict2tcl import dict2tcl

ROOT.gSystem.Load("libDelphes")

C_LIGHT = 2.99792458e8
C_LIGHT_MM_PER_NS = C_LIGHT * 1.0e-6

REPO_ROOT = Path(__file__).parent.parent
TESTS_DIR = Path(__file__).parent
DATA_DIR = TESTS_DIR / "data"
CARDS_DIR = REPO_ROOT / "cards"
MINBIAS_FILE = str(REPO_ROOT / "MinBias.pileup")

pythia8_reader_available = hasattr(ROOT, "DelphesPythia8Reader")
pythia8_pileup_available = hasattr(ROOT, "PileUpMergerPythia8")


@pytest.fixture(scope="function")
def load_delphes():
    modules = []
    refs = []

    def load(config, fout=0):
        conf_reader = ROOT.ExRootConfReader()
        if isinstance(config, dict):
            data = dict2tcl(config).encode()
            conf_reader.ReadData(".", data, len(data))
        else:
            conf_reader.ReadFile(config)
        writer = ROOT.ExRootTreeWriter(fout, "Delphes")
        module = ROOT.Delphes("Delphes")
        module.SetConfReader(conf_reader)
        module.SetTreeWriter(writer)
        factory = module.GetFactory()
        modules.append(module)
        refs.extend([conf_reader])
        return module, factory, writer

    yield load

    for m in modules:
        m.Finish()
    modules.clear()
    refs.clear()


def make_config(class_name, input_array="Delphes/inputParticles", output_array="outputParticles", **kwargs):
    return {
        "RandomSeed": 42,
        "ExecutionPath": ["TestModule"],
        "TestModule": {
            "Class": class_name,
            "InputArray": input_array,
            "OutputArray": output_array,
            **kwargs,
        },
    }


def build_config(class_name, defaults=None, **extra):
    params = dict(defaults or {})
    params.update(extra)
    return make_config(class_name, **params)


def make_candidate(factory, pt, eta, phi=0.0, energy=None, pid=0, charge=0, status=1, m1=-1):
    if energy is None:
        energy = pt * math.cosh(eta)
    c = factory.NewCandidate()
    c.Momentum.SetPtEtaPhiE(pt, eta, phi, energy)
    c.Position.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0)
    c.PID = pid
    c.Charge = charge
    c.Status = status
    c.M1 = m1
    return c


def make_jet(factory, pt, eta, phi=0.0, energy=None):
    return make_candidate(factory, pt, eta, phi, energy, pid=0, charge=0)


def make_parton(factory, pt, eta, pid, status=3, d1=-1, d2=-1, phi=0.0, energy=None):
    charge = 1 if pid > 0 else -1
    c = make_candidate(factory, pt, eta, phi, energy, pid=pid, charge=charge, status=status)
    c.D1 = d1
    c.D2 = d2
    return c


def make_particle(factory, pt, eta, pid, status=1, phi=0.0, energy=None):
    charge = 1 if pid > 0 else -1
    return make_candidate(factory, pt, eta, phi, energy, pid=pid, charge=charge, status=status)


def make_vertex_finder_track(factory, pt, eta, dz, error_dz, is_pu):
    c = make_candidate(factory, pt, eta, charge=1, pid=211)
    c.Position.SetXYZT(0.0, 0.0, dz * 1000.0, 0.0)
    c.InitialPosition.SetXYZT(0.0, 0.0, dz * 1000.0, 0.0)
    c.DZ = dz * 1000.0
    c.ErrorDZ = error_dz * 1000.0
    c.D0 = 0.0
    c.ErrorD0 = 0.001
    c.P = pt
    c.CtgTheta = 1.0
    c.Phi = 0.0
    c.IsPU = is_pu
    return c


def make_vertex(factory, x=0.0, y=0.0, z=0.0, t=0.0, is_pu=0):
    c = make_candidate(factory, 0.0, 0.0)
    c.Position.SetXYZT(x, y, z, t)
    c.IsPU = is_pu
    return c


def make_flavor_jet(factory, pt, eta, phi=0.0, flavor=0):
    jet = make_jet(factory, pt, eta, phi)
    jet.Flavor = flavor
    jet.FlavorAlgo = flavor
    jet.FlavorPhys = flavor
    return jet


CANDIDATE_DEFAULTS = {"pt": 50.0, "eta": 0.5, "phi": 0.0, "pid": 0, "charge": 0, "status": 1}
CANDIDATE_KEYS = ("pt", "eta", "phi", "pid", "charge", "status")

JET_DEFAULTS = {"pt": 50.0, "eta": 0.5, "phi": 0.0, "flavor": 0}
JET_KEYS = ("pt", "eta", "phi", "flavor")

PARTON_DEFAULTS = {"pt": 50.0, "eta": 0.5, "phi": 0.0, "pid": 0, "status": 3, "d1": -1, "d2": -1}
PARTON_KEYS = ("pt", "eta", "phi", "pid", "status", "d1", "d2")

PARTICLE_DEFAULTS = {"pt": 50.0, "eta": 0.5, "phi": 0.0, "pid": 0, "status": 1}
PARTICLE_KEYS = ("pt", "eta", "phi", "pid", "status")


def add_candidates(module, factory, array_name, specs, builder, defaults, consumed):
    array = module.ExportArray(array_name)
    for spec in specs:
        kwargs = dict(defaults)
        kwargs.update(spec)
        build_kwargs = {key: kwargs.pop(key) for key in consumed}
        candidate = builder(factory, **build_kwargs)
        for key, value in kwargs.items():
            setattr(candidate, key, value)
        array.Add(candidate)
    return array


def add_input_candidates(module, factory, specs):
    return add_candidates(module, factory, "inputParticles", specs, make_candidate, CANDIDATE_DEFAULTS, CANDIDATE_KEYS)


def add_generator_particles(module, factory, specs):
    return add_candidates(module, factory, "inputParticles", specs, make_particle, PARTICLE_DEFAULTS, PARTICLE_KEYS)


def add_input_jets(module, factory, specs):
    return add_candidates(module, factory, "inputJets", specs, make_flavor_jet, JET_DEFAULTS, JET_KEYS)


def add_input_partons(module, factory, specs):
    return add_candidates(module, factory, "inputPartons", specs, make_parton, PARTON_DEFAULTS, PARTON_KEYS)


VERTEX_FINDER_TRACKS = [
    (50.0, 0.5, 0.0, 0.001, 0),
    (40.0, 0.5, 0.001, 0.001, 0),
    (30.0, 0.5, -0.001, 0.001, 0),
    (20.0, 0.5, 0.002, 0.001, 0),
]


def run_vertex_finder_test(run_generic, config, tracks):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        for pt, eta, dz, error_dz, is_pu in tracks:
            input_array.Add(make_vertex_finder_track(factory, pt, eta, dz, error_dz, is_pu))

    return run_generic(config, setup=setup, outputs=("TestModule/tracks", "TestModule/vertices"))


@pytest.fixture(scope="function")
def run_generic(load_delphes):
    def run(config, setup=None, outputs=("TestModule/outputParticles",)):
        module, factory, _ = load_delphes(config)
        result = setup(module, factory) if setup else None
        module.Init()
        module.Process()
        if result is not None:
            return result
        if isinstance(outputs, dict):
            return {name: module.ImportArray(path) for name, path in outputs.items()}
        arrays = [module.ImportArray(name) for name in outputs]
        if len(arrays) == 1:
            return arrays[0]
        return arrays

    yield run


@pytest.fixture(scope="function")
def run_module(run_generic):
    def run(config, candidates):
        def setup(module, factory):
            add_input_candidates(module, factory, candidates)

        return run_generic(config, setup=setup)

    yield run


def run_repeated(load_delphes, config, n, setup, output="TestModule/outputParticles"):
    results = []
    for i in range(n):
        cfg = dict(config)
        cfg["RandomSeed"] = 42 + i
        module, factory, _ = load_delphes(cfg)
        setup(module, factory)
        module.Init()
        module.Process()
        results.append(module.ImportArray(output))
    return results


def mean(values):
    return sum(values) / len(values)


def stddev(values):
    m = mean(values)
    return (sum((v - m) ** 2 for v in values) / (len(values) - 1)) ** 0.5


@pytest.fixture(scope="function")
def run_tagging(run_generic):
    def run(config, jets, partons=None, particles=None, tracks=None):
        def setup(module, factory):
            jet_array = add_input_jets(module, factory, jets)
            if partons is not None:
                add_input_partons(module, factory, partons)
            if particles is not None:
                add_generator_particles(module, factory, particles)
            if tracks is not None:
                add_candidates(
                    module, factory, "inputTracks", tracks, make_candidate, CANDIDATE_DEFAULTS, CANDIDATE_KEYS
                )
            return jet_array

        return run_generic(config, setup=setup)

    yield run


@pytest.fixture(scope="function")
def run_calorimeter(run_generic):
    def run(config, particles, tracks=None):
        def setup(module, factory):
            particle_array = module.ExportArray("inputParticles")
            for pt, eta, phi, pid in particles:
                c = make_candidate(factory, pt, eta, phi, pid=pid)
                c.Position.SetPtEtaPhiE(1.0, eta, phi, 0)
                particle_array.Add(c)

            track_array = module.ExportArray("inputTracks")
            if tracks:
                for pt, eta, phi, charge, pid, track_resolution in tracks:
                    c = make_candidate(factory, pt, eta, phi, pid=pid, charge=charge)
                    c.Position.SetPtEtaPhiE(1.0, eta, phi, 0)
                    c.TrackResolution = track_resolution
                    track_array.Add(c)

        return run_generic(config, setup=setup, outputs=("TestModule/towers",))

    yield run


def candidate_snapshots(array, fields=("PID", "Momentum", "Position")):
    snaps = []
    for i in range(array.GetEntries()):
        c = array.At(i)
        values = []
        for field in fields:
            if field == "Momentum":
                mom = c.Momentum
                values.extend((mom.Px(), mom.Py(), mom.Pz(), mom.E()))
            elif field == "Position":
                pos = c.Position
                values.extend((pos.X(), pos.Y(), pos.Z(), pos.T()))
            else:
                values.append(getattr(c, field))
        snaps.append(tuple(values))
    return tuple(snaps)


def deterministic_snapshot(result, extract):
    if extract is not None:
        return extract(result)
    if isinstance(result, (tuple, list)):
        return tuple(candidate_snapshots(arr) for arr in result)
    return candidate_snapshots(result)


def flatten_snapshot(snap):
    if not isinstance(snap, (tuple, list)):
        return [snap]
    flat = []
    for value in snap:
        flat.extend(flatten_snapshot(value))
    return flat


def assert_deterministic(runner, extract=None, abs_tol=None):
    snap_a = deterministic_snapshot(runner(), extract)
    snap_b = deterministic_snapshot(runner(), extract)
    if abs_tol is None:
        assert snap_a == snap_b
    else:
        assert flatten_snapshot(snap_a) == pytest.approx(flatten_snapshot(snap_b), abs=abs_tol)


def run_seeds(load_delphes, config, seeds, setup, output="TestModule/outputParticles"):
    results = []
    for seed in seeds:
        cfg = dict(config)
        cfg["RandomSeed"] = seed
        module, factory, _ = load_delphes(cfg)
        setup(module, factory)
        module.Init()
        module.Process()
        results.append(module.ImportArray(output))
    return results


def repeated_values(load_delphes, config, n, setup, extract, output="TestModule/outputParticles"):
    results = run_repeated(load_delphes, config, n, setup, output=output)
    return [extract(r) for r in results]


def run_pileup(load_delphes, config, particles):
    module, factory, _ = load_delphes(config)
    input_array = module.ExportArray("inputParticles")
    for pt, charge in particles:
        c = make_candidate(factory, pt, 0.5, pid=211, charge=charge)
        c.Position.SetXYZT(0.0, 0.0, 0.0, 0.0)
        input_array.Add(c)
    module.Init()
    module.Process()
    return module.ImportArray("TestModule/stableParticles"), module.ImportArray("TestModule/vertices")


def read_tree_branch(load_delphes, tmp_path, build, config=None, index=1):
    if config is None:
        config = {}
    path = tmp_path / f"output_{index}.root"
    fout = ROOT.TFile(str(path), "RECREATE")
    module, _, writer = load_delphes(config, fout)
    build(module, writer)
    writer.Fill()
    writer.Write()
    module.Finish()
    fout.Close()
    fin = ROOT.TFile(str(path))
    reader = ROOT.ExRootTreeReader(fin.Get("Delphes"))
    return reader, fin


@contextlib.contextmanager
def reader_context(load_delphes, format, data_file, config=None, fout=0):
    if config is None:
        config = {}
    module, factory, writer = load_delphes(config, fout)
    arrays = [module.ExportArray(name) for name in ("allParticles", "stableParticles", "partons")]
    reader = getattr(ROOT, f"Delphes{format}Reader")()
    reader.OpenInputFile(str(data_file))
    try:
        yield module, factory, writer, arrays, reader
    finally:
        reader.CloseInputFile()


EVENT_FIELD_ATTRIBUTES = {
    "number": "Number",
    "process_id": "ProcessID",
    "mpi": "MPI",
    "weight": "Weight",
    "scale": "Scale",
    "alpha_qed": "AlphaQED",
    "alpha_qcd": "AlphaQCD",
    "id1": "ID1",
    "id2": "ID2",
    "x1": "X1",
    "x2": "X2",
    "scale_pdf": "ScalePDF",
    "pdf1": "PDF1",
    "pdf2": "PDF2",
    "cross_section": "CrossSection",
    "cross_section_error": "CrossSectionError",
}

EXACT_EVENT_FIELDS = ("number", "process_id", "mpi", "id1", "id2")
ABS_TOLERANCE_EVENT_FIELDS = ("pdf1", "pdf2")


def check_event_fields(event, expected):
    for key, attribute in EVENT_FIELD_ATTRIBUTES.items():
        if key not in expected:
            continue
        value = getattr(event, attribute)
        if key in EXACT_EVENT_FIELDS:
            assert value == expected[key], key
        elif key in ABS_TOLERANCE_EVENT_FIELDS:
            assert value == pytest.approx(expected[key], abs=1e-6), key
        else:
            assert value == pytest.approx(expected[key], rel=1e-5), key


def run_cli(args, timeout=120):
    return subprocess.run(
        [sys.executable, "-m", "delphes.delphes", *args],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        timeout=timeout,
        check=False,
    )
