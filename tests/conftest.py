import math
import pytest
import ROOT

from delphes.dict2tcl import dict2tcl

ROOT.gSystem.Load("libDelphes")


@pytest.fixture(scope="function")
def load_delphes():
    modules = []
    refs = []

    def load(config):
        conf_reader = ROOT.ExRootConfReader()
        if isinstance(config, dict):
            data = dict2tcl(config).encode()
            conf_reader.ReadData(".", data, len(data))
        else:
            conf_reader.ReadFile(config)
        writer = ROOT.ExRootTreeWriter()
        module = ROOT.Delphes("Delphes")
        module.SetConfReader(conf_reader)
        module.SetTreeWriter(writer)
        factory = module.GetFactory()
        modules.append(module)
        refs.extend([conf_reader, writer])
        return module, factory

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
        module, factory = load_delphes(config)
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
