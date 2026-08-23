import math

import pytest
import ROOT
from conftest import C_LIGHT, make_candidate

from delphes.dict2tcl import dict2tcl


def load_conf(config):
    conf_reader = ROOT.ExRootConfReader()
    data = dict2tcl(config).encode()
    conf_reader.ReadData(".", data, len(data))
    return conf_reader


def write_tree(tmp_path, config, array_names, fill, n_events=1):
    output_file = tmp_path / "delphes_tree.root"
    fout = ROOT.TFile(str(output_file), "RECREATE")
    writer = ROOT.ExRootTreeWriter(fout, "Delphes")
    module = ROOT.Delphes("Delphes")
    module.SetConfReader(load_conf(config))
    module.SetTreeWriter(writer)
    factory = module.GetFactory()

    for name in array_names:
        module.ExportArray(name)
    try:
        module.Init()
        for event in range(n_events):
            if event > 0:
                module.Clear()
                writer.Clear()
            fill(module, factory, event)
            module.Process()
            writer.Fill()
        writer.Write()
    finally:
        module.Finish()
        fout.Close()
    return output_file


def open_tree(path):
    fin = ROOT.TFile(str(path))
    tree = fin.Get("Delphes")
    reader = ROOT.ExRootTreeReader(tree)
    return reader, fin


def read_first_entry(reader, *names):
    branches = [reader.UseBranch(name) for name in names]
    reader.ReadEntry(0)
    return branches if len(branches) > 1 else branches[0]


def add_particle(factory, array, pt, eta=0.5, pid=211, charge=1, status=1):
    c = make_candidate(factory, pt, eta, pid=pid, charge=charge, status=status)
    c.Position.SetXYZT(0.0, 0.0, 0.0, 0.0)
    array.Add(c)
    return c


def add_track(factory, array, pt, eta=0.5, d0=0.01):
    c = make_candidate(factory, pt, eta, pid=211, charge=1)
    c.Position.SetXYZT(100.0, 0.0, 50.0, 0.0)
    c.D0 = d0
    c.AddCandidate(make_candidate(factory, 0.0, 0.0))
    array.Add(c)
    return c


def add_tower(factory, array, pt, eta=0.5, eem=10.0, ehad=20.0):
    c = make_candidate(factory, pt, eta, pid=22, charge=0)
    c.Eem = eem
    c.Ehad = ehad
    c.Edges[0] = eta - 0.05
    c.Edges[1] = eta + 0.05
    c.Edges[2] = -0.1
    c.Edges[3] = 0.1
    array.Add(c)


def add_jet(factory, array, pt, eta=0.5, flavor=5):
    c = make_candidate(factory, pt, eta)
    c.Flavor = flavor
    c.FlavorAlgo = flavor
    c.FlavorPhys = flavor
    array.Add(c)
    return c


def add_single(factory, array, value):
    c = make_candidate(factory, value, 0.0)
    array.Add(c)


def tree_writer_config(branches):
    return {
        "RandomSeed": 42,
        "ExecutionPath": ["TreeWriter"],
        "TreeWriter": {"Class": "TreeWriter", "Branch": branches},
    }


def test_roundtrip_particle_branch(tmp_path):
    config = tree_writer_config(["Delphes/inputParticles", "Particle", "GenParticle"])

    def fill(module, factory, event):
        array = module.ImportArray("Delphes/inputParticles")
        add_particle(factory, array, 50.0, pid=211, charge=1)
        add_particle(factory, array, 30.0, pid=11, charge=-1)

    path = write_tree(tmp_path, config, ["inputParticles"], fill)
    reader, fin = open_tree(path)
    assert reader.GetEntries() == 1
    branch = read_first_entry(reader, "Particle")
    assert branch.GetEntries() == 2
    p0 = branch.At(0)
    assert p0.PID == 211
    assert p0.Status == 1
    assert p0.Charge == 1
    assert p0.PT == pytest.approx(50.0, rel=1e-6)
    assert p0.Eta == pytest.approx(0.5, rel=1e-6)
    p1 = branch.At(1)
    assert p1.PID == 11
    assert p1.PT == pytest.approx(30.0, rel=1e-6)
    fin.Close()


def test_roundtrip_multiple_branch_types(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputParticles",
            "Particle",
            "GenParticle",
            "Delphes/inputTracks",
            "Track",
            "Track",
            "Delphes/inputTowers",
            "Tower",
            "Tower",
            "Delphes/inputJets",
            "Jet",
            "Jet",
            "Delphes/inputMissingET",
            "MissingET",
            "MissingET",
            "Delphes/inputRho",
            "Rho",
            "Rho",
            "Delphes/inputWeight",
            "Weight",
            "Weight",
        ]
    )

    array_names = (
        "inputParticles",
        "inputTracks",
        "inputTowers",
        "inputJets",
        "inputMissingET",
        "inputRho",
        "inputWeight",
    )

    def fill(module, factory, event):
        add_particle(factory, module.ImportArray("Delphes/inputParticles"), 50.0, pid=211, charge=1)
        add_track(factory, module.ImportArray("Delphes/inputTracks"), 40.0, d0=0.02)
        add_tower(factory, module.ImportArray("Delphes/inputTowers"), 25.0, eem=10.0, ehad=15.0)
        add_jet(factory, module.ImportArray("Delphes/inputJets"), 100.0, flavor=5)
        add_single(factory, module.ImportArray("Delphes/inputMissingET"), 35.0)
        add_single(factory, module.ImportArray("Delphes/inputRho"), 1.5)
        add_single(factory, module.ImportArray("Delphes/inputWeight"), 0.75)

    path = write_tree(tmp_path, config, array_names, fill)
    reader, fin = open_tree(path)
    particles, tracks, towers, jets, met, rho, weight = read_first_entry(
        reader, "Particle", "Track", "Tower", "Jet", "MissingET", "Rho", "Weight"
    )

    assert particles.GetEntries() == 1
    assert particles.At(0).PID == 211
    assert particles.At(0).PT == pytest.approx(50.0, rel=1e-6)

    assert tracks.GetEntries() == 1
    assert tracks.At(0).PT == pytest.approx(40.0, rel=1e-6)
    assert tracks.At(0).D0 == pytest.approx(0.02, rel=1e-6)

    assert towers.GetEntries() == 1
    assert towers.At(0).ET == pytest.approx(25.0, rel=1e-6)
    assert towers.At(0).Eem == pytest.approx(10.0, rel=1e-6)
    assert towers.At(0).Ehad == pytest.approx(15.0, rel=1e-6)
    assert towers.At(0).Edges[0] == pytest.approx(0.45, rel=1e-6)

    assert jets.GetEntries() == 1
    assert jets.At(0).PT == pytest.approx(100.0, rel=1e-6)
    assert jets.At(0).Flavor == 5

    assert met.GetEntries() == 1
    assert met.At(0).MET == pytest.approx(35.0, rel=1e-6)

    assert rho.GetEntries() == 1
    assert rho.At(0).Rho == pytest.approx(1.5, rel=1e-6)

    assert weight.GetEntries() == 1
    assert weight.At(0).Weight == pytest.approx(0.75, rel=1e-6)
    fin.Close()


def test_two_events(tmp_path):
    config = tree_writer_config(["Delphes/inputParticles", "Particle", "GenParticle"])

    def fill(module, factory, event):
        array = module.ImportArray("Delphes/inputParticles")
        add_particle(factory, array, 10.0 * (event + 1), pid=211, charge=1)

    path = write_tree(tmp_path, config, ["inputParticles"], fill, n_events=2)
    reader, fin = open_tree(path)
    assert reader.GetEntries() == 2
    branch = read_first_entry(reader, "Particle")
    assert branch.GetEntries() == 1
    assert branch.At(0).PT == pytest.approx(10.0, rel=1e-6)
    reader.ReadEntry(1)
    assert branch.GetEntries() == 1
    assert branch.At(0).PT == pytest.approx(20.0, rel=1e-6)
    fin.Close()


def test_empty_event(tmp_path):
    config = tree_writer_config(["Delphes/inputParticles", "Particle", "GenParticle"])

    def fill(module, factory, event):
        pass

    path = write_tree(tmp_path, config, ["inputParticles"], fill)
    reader, fin = open_tree(path)
    assert reader.GetEntries() == 1
    branch = read_first_entry(reader, "Particle")
    assert branch.GetEntries() == 0
    fin.Close()


def test_chained_execution_path(tmp_path):
    config = {
        "RandomSeed": 42,
        "ExecutionPath": ["EnergyScale", "Cloner", "TreeWriter"],
        "EnergyScale": {
            "Class": "EnergyScale",
            "InputArray": "Delphes/inputParticles",
            "OutputArray": "scaled",
            "ScaleFormula": 0.5,
        },
        "Cloner": {"Class": "Cloner", "InputArray": "EnergyScale/scaled", "OutputArray": "clones"},
        "TreeWriter": {
            "Class": "TreeWriter",
            "Branch": ["Cloner/clones", "Particle", "GenParticle"],
        },
    }

    def fill(module, factory, event):
        array = module.ImportArray("Delphes/inputParticles")
        add_particle(factory, array, 50.0, pid=211, charge=1)
        add_particle(factory, array, 30.0, pid=211, charge=1)

    path = write_tree(tmp_path, config, ["inputParticles"], fill)
    reader, fin = open_tree(path)
    branch = read_first_entry(reader, "Particle")
    assert branch.GetEntries() == 2

    assert branch.At(0).PT == pytest.approx(25.0, rel=1e-6)
    assert branch.At(1).PT == pytest.approx(15.0, rel=1e-6)
    fin.Close()


def test_unknown_branch_class_is_skipped(tmp_path, capfd):
    config = {
        "RandomSeed": 42,
        "ExecutionPath": ["TreeWriter"],
        "TreeWriter": {
            "Class": "TreeWriter",
            "Branch": [
                "Delphes/inputParticles",
                "Particle",
                "GenParticle",
                "Delphes/inputJets",
                "BadJets",
                "NoSuchClass",
                "Delphes/inputRho",
                "BadRho",
                "LHEFEvent",
            ],
        },
    }

    def fill(module, factory, event):
        add_particle(factory, module.ImportArray("Delphes/inputParticles"), 50.0, pid=211, charge=1)
        add_jet(factory, module.ImportArray("Delphes/inputJets"), 100.0, flavor=5)
        add_single(factory, module.ImportArray("Delphes/inputRho"), 1.5)

    path = write_tree(tmp_path, config, ["inputParticles", "inputJets", "inputRho"], fill)
    captured = capfd.readouterr()
    assert "** ERROR: cannot find class 'NoSuchClass'" in captured.out
    assert "** ERROR: cannot create branch for class 'LHEFEvent'" in captured.out

    reader, fin = open_tree(path)

    branch = read_first_entry(reader, "Particle")
    assert branch.GetEntries() == 1
    assert branch.At(0).PT == pytest.approx(50.0, rel=1e-6)

    tree = fin.Get("Delphes")
    branch_list = tree.GetListOfBranches()
    branch_names = [branch_list.At(i).GetName() for i in range(branch_list.GetEntries())]
    assert "Particle" in branch_names
    assert "BadJets" not in branch_names
    assert "BadRho" not in branch_names
    fin.Close()


def test_info_param_lands_in_tree_user_info(tmp_path, capfd):
    config = {
        "RandomSeed": 42,
        "ExecutionPath": ["TreeWriter"],
        "TreeWriter": {
            "Class": "TreeWriter",
            "Branch": ["Delphes/inputParticles", "Particle", "GenParticle"],
            "Info": ["TestInfo", 42.5],
        },
    }

    def fill(module, factory, event):
        add_particle(factory, module.ImportArray("Delphes/inputParticles"), 50.0, pid=211, charge=1)

    path = write_tree(tmp_path, config, ["inputParticles"], fill)
    fin = ROOT.TFile(str(path))
    tree = fin.Get("Delphes")
    user_info = tree.GetUserInfo()
    param = user_info.FindObject("TestInfo")
    assert param is not None

    assert param.IsA().GetName() == "TParameter<double>"

    ROOT.gROOT.ProcessLine(
        'printf("__INFOVAL__ %.17g\\n", '
        '((TParameter<double>*)(((TTree*)gDirectory->Get("Delphes"))->GetUserInfo()->'
        'FindObject("TestInfo")))->GetVal()); fflush(stdout);'
    )
    out = capfd.readouterr().out
    lines = [ln for ln in out.splitlines() if ln.startswith("__INFOVAL__")]
    assert len(lines) == 1
    assert float(lines[0].split()[1]) == pytest.approx(42.5, rel=1e-12)

    branch_list = tree.GetListOfBranches()
    branch_names = [branch_list.At(i).GetName() for i in range(branch_list.GetEntries())]
    assert "TestInfo" not in branch_names
    fin.Close()


def test_missing_et_fields_and_first_entry_only(tmp_path):
    config = tree_writer_config(["Delphes/inputMissingET", "MissingET", "MissingET"])

    def fill(module, factory, event):
        array = module.ImportArray("Delphes/inputMissingET")
        c1 = make_candidate(factory, 0.0, 0.0)
        c1.Momentum.SetPxPyPzE(3.0, 4.0, 0.0, 5.0)
        array.Add(c1)
        c2 = make_candidate(factory, 0.0, 0.0)
        c2.Momentum.SetPxPyPzE(0.0, 0.0, 7.0, 7.0)
        array.Add(c2)

    path = write_tree(tmp_path, config, ["inputMissingET"], fill)
    reader, fin = open_tree(path)
    branch = read_first_entry(reader, "MissingET")
    assert branch.GetEntries() == 1
    met = branch.At(0)
    assert met.MET == pytest.approx(5.0, rel=1e-6)
    assert met.Eta == pytest.approx(0.0, abs=1e-6)
    assert met.Phi == pytest.approx(math.atan2(-4.0, -3.0), rel=1e-5)
    fin.Close()


def test_scalar_ht_and_weight_first_entry_only(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputScalarHT",
            "ScalarHT",
            "ScalarHT",
            "Delphes/inputWeight",
            "Weight",
            "Weight",
        ]
    )

    def fill(module, factory, event):
        add_single(factory, module.ImportArray("Delphes/inputScalarHT"), 30.0)
        add_single(factory, module.ImportArray("Delphes/inputScalarHT"), 40.0)
        add_single(factory, module.ImportArray("Delphes/inputWeight"), 0.5)
        add_single(factory, module.ImportArray("Delphes/inputWeight"), 0.75)

    path = write_tree(tmp_path, config, ["inputScalarHT", "inputWeight"], fill)
    reader, fin = open_tree(path)
    ht_branch, weight_branch = read_first_entry(reader, "ScalarHT", "Weight")
    assert ht_branch.GetEntries() == 1
    assert ht_branch.At(0).HT == pytest.approx(30.0, rel=1e-6)
    assert weight_branch.GetEntries() == 1
    assert weight_branch.At(0).Weight == pytest.approx(0.5, rel=1e-6)
    fin.Close()


def test_vertex_branch_roundtrip(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputParticles",
            "Particle",
            "GenParticle",
            "Delphes/inputVertices",
            "Vertex",
            "Vertex",
        ]
    )

    def fill(module, factory, event):
        gen1 = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 10.0, pid=211, charge=1)
        gen2 = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 20.0, pid=22, charge=0)
        vertex = make_candidate(factory, 0.0, 0.0)
        vertex.Position.SetXYZT(1.0, 2.0, 3.0, 4000.0)
        vertex.PositionError.SetXYZT(0.1, 0.2, 0.3, 10.0)
        vertex.ClusterIndex = 7
        vertex.ClusterNDF = 12
        vertex.ClusterSigma = 0.5
        vertex.SumPT2 = 100.0
        vertex.BTVSumPT2 = 50.0
        vertex.GenDeltaZ = 1.5
        vertex.GenSumPT2 = 90.0
        vertex.AddCandidate(gen1)
        vertex.AddCandidate(gen2)
        module.ImportArray("Delphes/inputVertices").Add(vertex)

    path = write_tree(tmp_path, config, ["inputParticles", "inputVertices"], fill)
    reader, fin = open_tree(path)

    reader.UseBranch("Particle")
    branch = read_first_entry(reader, "Vertex")
    assert branch.GetEntries() == 1
    v = branch.At(0)
    assert v.Index == 7
    assert v.NDF == 12
    assert v.Sigma == pytest.approx(0.5, rel=1e-6)
    assert v.SumPT2 == pytest.approx(100.0, rel=1e-6)
    assert v.BTVSumPT2 == pytest.approx(50.0, rel=1e-6)
    assert v.GenDeltaZ == pytest.approx(1.5, rel=1e-6)
    assert v.GenSumPT2 == pytest.approx(90.0, rel=1e-6)
    assert v.X == pytest.approx(1.0, rel=1e-6)
    assert v.Y == pytest.approx(2.0, rel=1e-6)
    assert v.Z == pytest.approx(3.0, rel=1e-6)
    assert v.T == pytest.approx(4.0 / C_LIGHT, rel=1e-5)
    assert v.ErrorX == pytest.approx(0.1, rel=1e-6)
    assert v.ErrorY == pytest.approx(0.2, rel=1e-6)
    assert v.ErrorZ == pytest.approx(0.3, rel=1e-6)
    assert v.ErrorT == pytest.approx(10.0e-3 / C_LIGHT, rel=1e-5)

    assert v.Constituents.GetEntries() == 2
    pids = sorted(v.Constituents.At(i).PID for i in range(2))
    pts = sorted(v.Constituents.At(i).PT for i in range(2))
    assert pids == [22, 211]
    assert pts == [pytest.approx(10.0, rel=1e-6), pytest.approx(20.0, rel=1e-6)]
    fin.Close()


def test_jet_constituents_refarray_readback(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputParticles",
            "Particle",
            "GenParticle",
            "Delphes/inputJets",
            "Jet",
            "Jet",
        ]
    )

    def fill(module, factory, event):
        gen1 = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 10.0, pid=211, charge=1)
        gen2 = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 20.0, pid=22, charge=0)
        jet = add_jet(factory, module.ImportArray("Delphes/inputJets"), 30.0, flavor=5)
        jet.AddCandidate(gen1)
        jet.AddCandidate(gen2)

    path = write_tree(tmp_path, config, ["inputParticles", "inputJets"], fill)
    reader, fin = open_tree(path)

    reader.UseBranch("Particle")
    branch = read_first_entry(reader, "Jet")
    assert branch.GetEntries() == 1
    jet = branch.At(0)
    assert jet.PT == pytest.approx(30.0, rel=1e-6)
    assert jet.Eta == pytest.approx(0.5, rel=1e-6)
    assert jet.Mass == pytest.approx(0.0, abs=1e-6)
    assert jet.Flavor == 5
    assert jet.EhadOverEem == pytest.approx(999.9, rel=1e-6)
    assert jet.Constituents.GetEntries() == 2
    pids = sorted(jet.Constituents.At(i).PID for i in range(2))
    pts = sorted(jet.Constituents.At(i).PT for i in range(2))
    assert pids == [22, 211]
    assert pts == [pytest.approx(10.0, rel=1e-6), pytest.approx(20.0, rel=1e-6)]

    assert jet.Particles.GetEntries() == 2
    fin.Close()


def test_track_particle_ref_readback(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputParticles",
            "Particle",
            "GenParticle",
            "Delphes/inputTracks",
            "Track",
            "Track",
        ]
    )

    def fill(module, factory, event):
        gen = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 25.0, pid=211, charge=1)
        track = make_candidate(factory, 25.0, 0.5, pid=211, charge=1)
        track.Position.SetXYZT(100.0, 0.0, 50.0, 0.0)
        track.InitialPosition.SetXYZT(1.0, 2.0, 3.0, 4000.0)
        track.D0 = 0.02
        track.ErrorD0 = 0.001
        track.DZ = -0.05
        track.ErrorDZ = 0.002
        track.ErrorP = 0.1
        track.ErrorPT = 0.15
        track.ErrorC = 0.0005
        track.C = 0.01
        track.L = 1.2
        track.Nclusters = 10.0
        track.dNdx = 1500.0
        track.Xd = 0.5
        track.Yd = 0.6
        track.Zd = 0.7
        track.XFirstHit = 1.1
        track.YFirstHit = 1.2
        track.ZFirstHit = 1.3
        track.ErrorT = 1.0
        track.IsPU = 1
        track.ClusterIndex = 3
        track.TrackCovariance[0][1] = 1e-7
        track.TrackCovariance[2][3] = 5e-5
        track.TrackCovariance[0][4] = 2e-8
        track.AddCandidate(gen)
        module.ImportArray("Delphes/inputTracks").Add(track)

    path = write_tree(tmp_path, config, ["inputParticles", "inputTracks"], fill)
    reader, fin = open_tree(path)

    reader.UseBranch("Particle")
    branch = read_first_entry(reader, "Track")
    assert branch.GetEntries() == 1
    t = branch.At(0)

    assert t.PT == pytest.approx(25.0, rel=1e-6)
    assert t.P == pytest.approx(25.0 * math.cosh(0.5), rel=1e-6)
    assert t.Eta == pytest.approx(0.5, rel=1e-6)
    assert t.Phi == pytest.approx(0.0, abs=1e-6)
    assert t.CtgTheta == pytest.approx(math.sinh(0.5), rel=1e-6)
    assert t.C == pytest.approx(0.01, rel=1e-6)
    assert t.Mass == pytest.approx(0.0, abs=1e-6)

    assert t.XOuter == pytest.approx(100.0, rel=1e-6)
    assert t.ZOuter == pytest.approx(50.0, rel=1e-6)
    assert t.EtaOuter == pytest.approx(math.atanh(50.0 / math.sqrt(12500.0)), rel=1e-6)
    assert t.X == pytest.approx(1.0, rel=1e-6)
    assert t.Y == pytest.approx(2.0, rel=1e-6)
    assert t.Z == pytest.approx(3.0, rel=1e-6)
    assert t.T == pytest.approx(4.0 / C_LIGHT, rel=1e-5)
    assert t.ErrorT == pytest.approx(1.0e-3 / C_LIGHT, rel=1e-5)

    assert t.D0 == pytest.approx(0.02, rel=1e-6)
    assert t.ErrorD0 == pytest.approx(0.001, rel=1e-6)
    assert t.DZ == pytest.approx(-0.05, rel=1e-6)
    assert t.ErrorDZ == pytest.approx(0.002, rel=1e-6)
    assert t.ErrorP == pytest.approx(0.1, rel=1e-6)
    assert t.ErrorPT == pytest.approx(0.15, rel=1e-6)
    assert t.L == pytest.approx(1.2, rel=1e-6)
    assert t.Nclusters == pytest.approx(10.0, rel=1e-6)
    assert t.dNdx == pytest.approx(1500.0, rel=1e-6)
    assert t.Xd == pytest.approx(0.5, rel=1e-6)
    assert t.ZFirstHit == pytest.approx(1.3, rel=1e-6)

    assert t.ErrorD0Phi == pytest.approx(1e-7 * 1e3, rel=1e-6)
    assert t.ErrorCDZ == pytest.approx(5e-5, rel=1e-6)
    assert t.ErrorD0CtgTheta == pytest.approx(2e-8 * 1e3, rel=1e-6)

    assert t.IsPU == 1
    assert t.HardEnergyFraction == pytest.approx(0.0, abs=1e-6)
    assert t.VertexIndex == 3

    particle = t.Particle.GetObject()
    assert particle is not None
    assert particle.PID == 211
    assert particle.PT == pytest.approx(25.0, rel=1e-6)
    fin.Close()


def test_photon_branch_roundtrip(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputParticles",
            "Particle",
            "GenParticle",
            "Delphes/inputPhotons",
            "Photon",
            "Photon",
        ]
    )

    def fill(module, factory, event):
        gen1 = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 10.0, pid=211, charge=1)
        gen2 = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 20.0, pid=22, charge=0)
        array = module.ImportArray("Delphes/inputPhotons")
        p1 = make_candidate(factory, 30.0, 0.5, pid=22, charge=0)
        p1.Eem = 4.0
        p1.Ehad = 6.0
        p1.Position.SetXYZT(5.0, 6.0, 7.0, 1000.0)
        p1.PositionError.SetXYZT(0.01, 0.02, 0.03, 0.04)
        p1.DecayPosition.SetXYZT(1.5, 2.5, 3.5, 0.0)
        p1.IsolationVar = 0.15
        p1.IsolationVarRhoCorr = 0.1
        p1.SumPtCharged = 0.5
        p1.SumPtNeutral = 0.3
        p1.SumPtChargedPU = 0.05
        p1.SumPt = 0.8
        p1.Status = 1
        p1.AddCandidate(gen1)
        array.Add(p1)
        p2 = make_candidate(factory, 20.0, -0.3, pid=22, charge=0)
        p2.Eem = 0.0
        p2.Ehad = 3.0
        p2.Status = 2
        p2.AddCandidate(gen2)
        array.Add(p2)

    path = write_tree(tmp_path, config, ["inputParticles", "inputPhotons"], fill)
    reader, fin = open_tree(path)

    reader.UseBranch("Particle")
    branch = read_first_entry(reader, "Photon")
    assert branch.GetEntries() == 2

    p1, p2 = branch.At(0), branch.At(1)
    assert p1.PT == pytest.approx(30.0, rel=1e-6)
    assert p2.PT == pytest.approx(20.0, rel=1e-6)
    assert p1.Eta == pytest.approx(0.5, rel=1e-6)
    assert p1.E == pytest.approx(30.0 * math.cosh(0.5), rel=1e-6)
    assert p1.T == pytest.approx(1.0 / C_LIGHT, rel=1e-5)
    assert p1.X == pytest.approx(5.0, rel=1e-6)
    assert p1.Y == pytest.approx(6.0, rel=1e-6)
    assert p1.Z == pytest.approx(7.0, rel=1e-6)
    assert p1.ErrorTheta == pytest.approx(0.01, rel=1e-6)
    assert p1.ErrorPhi == pytest.approx(0.02, rel=1e-6)
    assert p1.ThetaP == pytest.approx(1.5, rel=1e-6)
    assert p1.PhiP == pytest.approx(2.5, rel=1e-6)
    assert p1.ErrorThetaP == pytest.approx(0.03, rel=1e-6)
    assert p1.ErrorPhiP == pytest.approx(0.04, rel=1e-6)
    assert p1.ErrorE == pytest.approx(3.5, rel=1e-6)
    assert p1.IsolationVar == pytest.approx(0.15, rel=1e-6)
    assert p1.IsolationVarRhoCorr == pytest.approx(0.1, rel=1e-6)
    assert p1.SumPtCharged == pytest.approx(0.5, rel=1e-6)
    assert p1.SumPtNeutral == pytest.approx(0.3, rel=1e-6)
    assert p1.SumPtChargedPU == pytest.approx(0.05, rel=1e-6)
    assert p1.SumPt == pytest.approx(0.8, rel=1e-6)
    assert p1.Status == 1
    assert p1.EhadOverEem == pytest.approx(6.0 / 4.0, rel=1e-6)
    assert p2.EhadOverEem == pytest.approx(999.9, rel=1e-6)
    assert p2.Status == 2

    assert p1.Particles.GetEntries() == 1
    assert p1.Particles.At(0).PT == pytest.approx(10.0, rel=1e-6)
    assert p2.Particles.GetEntries() == 1
    assert p2.Particles.At(0).PT == pytest.approx(20.0, rel=1e-6)
    fin.Close()


def test_electron_branch_roundtrip(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputParticles",
            "Particle",
            "GenParticle",
            "Delphes/inputElectrons",
            "Electron",
            "Electron",
        ]
    )

    def fill(module, factory, event):
        gen = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 15.0, pid=11, charge=-1)
        e = make_candidate(factory, 40.0, 1.0, pid=11, charge=-1)
        e.Eem = 5.0
        e.Ehad = 10.0
        e.D0 = 0.03
        e.ErrorD0 = 0.002
        e.DZ = 0.07
        e.ErrorDZ = 0.003
        e.IsolationVar = 0.2
        e.AddCandidate(gen)
        module.ImportArray("Delphes/inputElectrons").Add(e)

    path = write_tree(tmp_path, config, ["inputParticles", "inputElectrons"], fill)
    reader, fin = open_tree(path)

    reader.UseBranch("Particle")
    branch = read_first_entry(reader, "Electron")
    assert branch.GetEntries() == 1
    e = branch.At(0)
    assert e.PT == pytest.approx(40.0, rel=1e-6)
    assert e.Eta == pytest.approx(1.0, rel=1e-6)
    assert e.Charge == -1
    assert e.D0 == pytest.approx(0.03, rel=1e-6)
    assert e.ErrorD0 == pytest.approx(0.002, rel=1e-6)
    assert e.DZ == pytest.approx(0.07, rel=1e-6)
    assert e.ErrorDZ == pytest.approx(0.003, rel=1e-6)
    assert e.IsolationVar == pytest.approx(0.2, rel=1e-6)

    assert e.EhadOverEem == pytest.approx(0.0, abs=1e-6)
    particle = e.Particle.GetObject()
    assert particle is not None
    assert particle.PID == 11
    assert particle.PT == pytest.approx(15.0, rel=1e-6)
    fin.Close()


def test_muon_branch_roundtrip(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputParticles",
            "Particle",
            "GenParticle",
            "Delphes/inputMuons",
            "Muon",
            "Muon",
        ]
    )

    def fill(module, factory, event):
        gen = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 20.0, pid=13, charge=1)
        m = make_candidate(factory, 60.0, -0.8, pid=13, charge=1)
        m.D0 = 0.04
        m.ErrorD0 = 0.004
        m.DZ = -0.02
        m.ErrorDZ = 0.005
        m.SumPt = 1.25
        m.AddCandidate(gen)
        module.ImportArray("Delphes/inputMuons").Add(m)

    path = write_tree(tmp_path, config, ["inputParticles", "inputMuons"], fill)
    reader, fin = open_tree(path)

    reader.UseBranch("Particle")
    branch = read_first_entry(reader, "Muon")
    assert branch.GetEntries() == 1
    m = branch.At(0)
    assert m.PT == pytest.approx(60.0, rel=1e-6)
    assert m.Eta == pytest.approx(-0.8, rel=1e-6)
    assert m.Charge == 1
    assert m.D0 == pytest.approx(0.04, rel=1e-6)
    assert m.ErrorD0 == pytest.approx(0.004, rel=1e-6)
    assert m.DZ == pytest.approx(-0.02, rel=1e-6)
    assert m.ErrorDZ == pytest.approx(0.005, rel=1e-6)
    assert m.SumPt == pytest.approx(1.25, rel=1e-6)
    particle = m.Particle.GetObject()
    assert particle is not None
    assert particle.PID == 13
    assert particle.PT == pytest.approx(20.0, rel=1e-6)
    fin.Close()


def test_particle_flow_candidate_branch_roundtrip(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputParticles",
            "Particle",
            "GenParticle",
            "Delphes/inputPFC",
            "ParticleFlowCandidate",
            "ParticleFlowCandidate",
        ]
    )

    def fill(module, factory, event):
        gen = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 8.0, pid=221)
        array = module.ImportArray("Delphes/inputPFC")
        pfc_c = make_candidate(factory, 30.0, 0.4, pid=211, charge=1)
        pfc_c.IsPU = 1
        pfc_c.D0 = 0.01
        pfc_c.DZ = 0.02
        pfc_c.Nclusters = 5.0
        pfc_c.dNdx = 1200.0
        pfc_c.L = 0.8
        pfc_c.ErrorP = 0.05
        pfc_c.Eem = 2.0
        pfc_c.Ehad = 3.0
        pfc_c.Etrk = 1.5
        pfc_c.Edges[0] = 0.35
        pfc_c.Edges[1] = 0.45
        pfc_c.Edges[2] = -0.05
        pfc_c.Edges[3] = 0.05
        pfc_c.NTimeHits = 3
        pfc_c.Position.SetXYZT(50.0, 10.0, 20.0, 500.0)
        pfc_c.AddCandidate(gen)
        array.Add(pfc_c)
        pfc_n = make_candidate(factory, 20.0, -0.6, pid=22, charge=0)
        pfc_n.BetaStar = 0.25
        array.Add(pfc_n)

    path = write_tree(tmp_path, config, ["inputParticles", "inputPFC"], fill)
    reader, fin = open_tree(path)

    reader.UseBranch("Particle")
    branch = read_first_entry(reader, "ParticleFlowCandidate")
    assert branch.GetEntries() == 2

    c, n = branch.At(0), branch.At(1)
    assert c.PT == pytest.approx(30.0, rel=1e-6)
    assert c.E == pytest.approx(30.0 * math.cosh(0.4), rel=1e-6)
    assert c.P == pytest.approx(30.0 * math.cosh(0.4), rel=1e-6)
    assert c.CtgTheta == pytest.approx(math.sinh(0.4), rel=1e-6)
    assert c.D0 == pytest.approx(0.01, rel=1e-6)
    assert c.DZ == pytest.approx(0.02, rel=1e-6)
    assert c.Nclusters == pytest.approx(5.0, rel=1e-6)
    assert c.dNdx == pytest.approx(1200.0, rel=1e-6)
    assert c.L == pytest.approx(0.8, rel=1e-6)
    assert c.ErrorP == pytest.approx(0.05, rel=1e-6)
    assert c.Eem == pytest.approx(2.0, rel=1e-6)
    assert c.Ehad == pytest.approx(3.0, rel=1e-6)
    assert c.Etrk == pytest.approx(1.5, rel=1e-6)
    assert c.Edges[0] == pytest.approx(0.35, rel=1e-6)
    assert c.Edges[3] == pytest.approx(0.05, rel=1e-6)
    assert c.NTimeHits == 3
    assert c.XOuter == pytest.approx(50.0, rel=1e-6)
    assert c.TOuter == pytest.approx(0.5 / C_LIGHT, rel=1e-5)

    assert c.IsPU == 1
    assert c.HardEnergyFraction == pytest.approx(0.0, abs=1e-6)

    assert n.PT == pytest.approx(20.0, rel=1e-6)
    assert n.Charge == 0
    assert n.HardEnergyFraction == pytest.approx(0.25, rel=1e-6)

    assert c.Particles.GetEntries() == 1
    assert c.Particles.At(0).PT == pytest.approx(8.0, rel=1e-6)
    assert n.Particles.GetEntries() == 0
    fin.Close()


def test_csc_cluster_branch_roundtrip(tmp_path):
    config = tree_writer_config(["Delphes/inputCscCluster", "CscCluster", "CscCluster"])

    def fill(module, factory, event):
        csc = make_candidate(factory, 0.0, 0.0)

        csc.Momentum.SetPxPyPzE(100.0, 0.0, 0.0, 150.0)
        csc.PID = 12
        csc.Eem = 10.0
        csc.Ehad = 20.0
        csc.DecayPosition.SetXYZT(30.0, 40.0, 50.0, 0.0)
        module.ImportArray("Delphes/inputCscCluster").Add(csc)

    path = write_tree(tmp_path, config, ["inputCscCluster"], fill)
    reader, fin = open_tree(path)
    branch = read_first_entry(reader, "CscCluster")
    assert branch.GetEntries() == 1
    c = branch.At(0)
    assert c.PT == pytest.approx(100.0, rel=1e-6)
    assert c.Px == pytest.approx(100.0, rel=1e-6)
    assert c.Pz == pytest.approx(0.0, abs=1e-6)
    assert c.E == pytest.approx(150.0, rel=1e-6)
    assert c.pid == 12
    assert c.Eem == pytest.approx(10.0, rel=1e-6)
    assert c.Ehad == pytest.approx(20.0, rel=1e-6)
    decay_distance = math.sqrt(30.0**2 + 40.0**2 + 50.0**2)
    beta = 100.0 / 150.0
    gamma = 1.0 / math.sqrt(1.0 - beta * beta)
    assert c.beta == pytest.approx(beta, rel=1e-5)
    assert c.ctau == pytest.approx(decay_distance / (beta * gamma), rel=1e-5)
    assert c.T == pytest.approx(decay_distance * (1.0 / beta - 1.0) * 1e-3 / C_LIGHT * 1e9, rel=1e-5)
    assert c.X == pytest.approx(30.0, rel=1e-6)
    assert c.Y == pytest.approx(40.0, rel=1e-6)
    assert c.Z == pytest.approx(50.0, rel=1e-6)
    assert c.R == pytest.approx(50.0, rel=1e-6)
    fin.Close()


def test_hector_hit_branch_roundtrip(tmp_path):
    config = tree_writer_config(
        [
            "Delphes/inputParticles",
            "Particle",
            "GenParticle",
            "Delphes/inputHectorHits",
            "HectorHit",
            "HectorHit",
        ]
    )

    def fill(module, factory, event):
        gen = add_particle(factory, module.ImportArray("Delphes/inputParticles"), 50.0, pid=2212, charge=1)
        hit = make_candidate(factory, 0.0, 0.0)
        hit.Momentum.SetPxPyPzE(0.03, -0.05, 100.0, 100.00001)
        hit.Position.SetXYZT(10.0, -20.0, 25.0, 80.0)
        hit.AddCandidate(gen)
        module.ImportArray("Delphes/inputHectorHits").Add(hit)

    path = write_tree(tmp_path, config, ["inputParticles", "inputHectorHits"], fill)
    reader, fin = open_tree(path)

    reader.UseBranch("Particle")
    branch = read_first_entry(reader, "HectorHit")
    assert branch.GetEntries() == 1
    h = branch.At(0)
    assert h.E == pytest.approx(100.00001, rel=1e-6)
    assert h.Tx == pytest.approx(0.03, rel=1e-6)
    assert h.Ty == pytest.approx(-0.05, rel=1e-6)

    assert h.T == pytest.approx(80.0, rel=1e-6)
    assert h.X == pytest.approx(10.0, rel=1e-6)
    assert h.Y == pytest.approx(-20.0, rel=1e-6)
    assert h.S == pytest.approx(25.0, rel=1e-6)
    particle = h.Particle.GetObject()
    assert particle is not None
    assert particle.PID == 2212
    assert particle.PT == pytest.approx(50.0, rel=1e-6)
    fin.Close()
