import struct

import pytest
import ROOT
from conftest import MINBIAS_FILE, make_config

from delphes.dict2tcl import dict2tcl


def expect_runtime_error(fn, message_part=None):
    with pytest.raises(ROOT.std.runtime_error) as exc_info:
        fn()
    if message_part is not None:
        assert message_part in str(exc_info.value)


def xdr_u32(value):
    return struct.pack(">I", value & 0xFFFFFFFF)


def xdr_u64(value):
    return struct.pack(">Q", value & 0xFFFFFFFFFFFFFFFF)


def xdr_string(text):
    raw = text.encode()
    padding = (4 - len(raw) % 4) % 4
    return xdr_u32(len(raw)) + raw + b"\x00" * padding


def write_stdhep(tmp_path, name, blob):
    path = tmp_path / name
    path.write_bytes(blob)
    return str(path)


def read_stdhep(factory, path):
    reader = ROOT.DelphesSTDHEPReader()
    reader.OpenInputFile(path)
    reader.ReadEvent(factory, ROOT.TObjArray(), ROOT.TObjArray(), ROOT.TObjArray())
    reader.CloseInputFile()


def test_init_without_conf_reader():
    module = ROOT.Delphes("Delphes")
    expect_runtime_error(module.Init, "can't access configuration reader")
    module.Finish()


def test_init_without_tree_writer():
    conf_reader = ROOT.ExRootConfReader()
    data = dict2tcl({"ExecutionPath": []}).encode()
    conf_reader.ReadData(".", data, len(data))
    module = ROOT.Delphes("Delphes")
    module.SetConfReader(conf_reader)
    expect_runtime_error(module.Init, "can't access tree writer")
    module.Finish()


def test_unconfigured_module_in_execution_path(load_delphes):
    module, _, _ = load_delphes({"ExecutionPath": ["MissingModule"]})
    expect_runtime_error(module.Init, "specified in ExecutionPath but not configured")


def test_unknown_module_class(load_delphes):
    module, _, _ = load_delphes(make_config("NoSuchModuleClass"))
    expect_runtime_error(module.Init, "can't find class 'NoSuchModuleClass'")


def test_module_class_not_inheriting_delphes_module(load_delphes):
    module, _, _ = load_delphes(make_config("TLorentzVector"))
    expect_runtime_error(module.Init, "does not inherit from DelphesModule")


def test_module_int_string(load_delphes):
    module, _, _ = load_delphes(make_config("ExampleModule", IntParam="abc"))
    expect_runtime_error(module.Init, "is not an integer")


def test_module_int_out_of_range(load_delphes):
    module, _, _ = load_delphes(make_config("ExampleModule", IntParam=2**40))
    expect_runtime_error(module.Init, "is not an integer")


def test_module_list_index_out_of_range(load_delphes):
    module, _, _ = load_delphes({"Delphes::ListParam": [0, 1, 2]})
    expect_runtime_error(lambda: module.GetInt("ListParam", 0, 3), "is out of range")


def test_module_missing_required_key(load_delphes):
    module, _, _ = load_delphes(make_config("ExampleModule", InputArray="Delphes/nonexistentArray"))
    expect_runtime_error(module.Init, "can't access input list 'Delphes/nonexistentArray'")


def test_invalid_efficiency_formula(load_delphes):
    module, _, _ = load_delphes(make_config("Efficiency", EfficiencyFormula="pt +"))
    expect_runtime_error(module.Init, "Invalid formula")


def test_invalid_conversion_map_formula(load_delphes):
    module, _, _ = load_delphes(make_config("PhotonConversions", ConversionMap="rho +"))
    expect_runtime_error(module.Init, "Invalid formula")


def test_invalid_csc_cluster_formula(load_delphes):
    module, _, _ = load_delphes(make_config("CscClusterEfficiency", EfficiencyFormula="decayR +"))
    expect_runtime_error(module.Init, "Invalid formula")


def test_invalid_vertex_distribution_formula(load_delphes):
    module, _, _ = load_delphes(
        make_config(
            "PileUpMerger",
            VertexDistributionFormula="z +",
            PileUpFile=MINBIAS_FILE,
            MeanPileUp=1,
            PileUpDistribution=2,
            InputArray="Delphes/inputParticles",
        )
    )
    module.ExportArray("inputParticles")
    expect_runtime_error(module.Init, "Invalid formula")


def test_jet_fake_particle_invalid_pdg(load_delphes):
    module, _, _ = load_delphes(make_config("JetFakeParticle", EfficiencyFormula=[999, "1.0"]))
    expect_runtime_error(module.Init, "Jets can only fake into electrons, muons or photons")


def test_weighter_too_many_pdg_codes(load_delphes):
    module, _, _ = load_delphes(make_config("Weighter", OutputArray="weight", Weight=[[1, 2, 3, 4, 5], 1.0]))
    expect_runtime_error(module.Init, "only 1, 2, 3 or 4 PDG codes")


def tau_out_of_range_config(module_class):
    return make_config(
        module_class,
        ParticleInputArray="Delphes/inputParticles",
        PartonInputArray="Delphes/partons",
        TrackInputArray="Delphes/tracks",
        JetInputArray="Delphes/jets",
    )


def test_tau_tagging_daughter_index_out_of_range(load_delphes):
    module, factory, _ = load_delphes(tau_out_of_range_config("TauTagging"))
    particles = module.ExportArray("inputParticles")
    p = factory.NewCandidate()
    p.Momentum.SetPtEtaPhiE(10.0, 0.0, 0.0, 10.0)
    particles.Add(p)

    partons = module.ExportArray("partons")
    tau = factory.NewCandidate()
    tau.Momentum.SetPtEtaPhiE(50.0, 0.5, 0.0, 60.0)
    tau.PID = 15
    tau.D1 = 2
    tau.D2 = 2
    partons.Add(tau)
    module.ExportArray("jets")

    module.Init()
    expect_runtime_error(module.Process, "tau's daughter index is greater than the ParticleInputArray size")


def test_track_counting_tau_tagging_daughter_index_out_of_range(load_delphes):
    module, factory, _ = load_delphes(tau_out_of_range_config("TrackCountingTauTagging"))
    particles = module.ExportArray("inputParticles")
    p = factory.NewCandidate()
    p.Momentum.SetPtEtaPhiE(10.0, 0.0, 0.0, 10.0)
    particles.Add(p)

    partons = module.ExportArray("partons")
    tau = factory.NewCandidate()
    tau.Momentum.SetPtEtaPhiE(50.0, 0.5, 0.0, 60.0)
    tau.PID = 15
    tau.D1 = 2
    tau.D2 = 2
    partons.Add(tau)
    module.ExportArray("tracks")
    module.ExportArray("jets")

    module.Init()
    expect_runtime_error(module.Process, "tau's daughter index is greater than the ParticleInputArray size")


def test_hepmc2_missing_input_file(tmp_path):
    reader = ROOT.DelphesHepMC2Reader()
    expect_runtime_error(lambda: reader.OpenInputFile(str(tmp_path / "missing.hepmc2")))
    reader.CloseInputFile()


def test_hepmc3_missing_input_file(tmp_path):
    reader = ROOT.DelphesHepMC3Reader()
    expect_runtime_error(lambda: reader.OpenInputFile(str(tmp_path / "missing.hepmc3")))
    reader.CloseInputFile()


def test_lhef_missing_input_file(tmp_path):
    reader = ROOT.DelphesLHEFReader()
    expect_runtime_error(lambda: reader.OpenInputFile(str(tmp_path / "missing.lhe")))
    reader.CloseInputFile()


def test_pythia8_bad_configuration_file(tmp_path):
    reader = ROOT.DelphesPythia8Reader()
    expect_runtime_error(
        lambda: reader.OpenInputFile(str(tmp_path / "missing.cmnd")),
        "can't read Pythia8 configuration file",
    )
    reader.CloseInputFile()


def test_stdhep_missing_input_file(tmp_path):
    reader = ROOT.DelphesSTDHEPReader()
    expect_runtime_error(lambda: reader.OpenInputFile(str(tmp_path / "missing.stdhep")))
    reader.CloseInputFile()


def test_stdhep_unsupported_block_type(load_delphes, tmp_path):
    _, factory, _ = load_delphes({})
    path = write_stdhep(tmp_path, "bad_block.stdhep", xdr_u32(999) + xdr_u32(0))
    expect_runtime_error(lambda: read_stdhep(factory, path), "Unsupported block type")


def test_stdhep_unknown_file_version(load_delphes, tmp_path):
    _, factory, _ = load_delphes({})
    blob = xdr_u32(1) + xdr_u32(0) + xdr_string("X")
    path = write_stdhep(tmp_path, "bad_version.stdhep", blob)
    expect_runtime_error(lambda: read_stdhep(factory, path), "Unknown file format version")


def test_stdhep_ntuples_not_supported(load_delphes, tmp_path):
    _, factory, _ = load_delphes({})
    blob = (
        xdr_u32(1)
        + xdr_u32(0)
        + xdr_string("2.00")
        + xdr_u32(0) * 3
        + xdr_u32(0)
        + xdr_u32(0)
        + b"\x00" * 8
        + xdr_u32(0)
        + xdr_u32(5)
    )
    path = write_stdhep(tmp_path, "ntuples.stdhep", blob)
    expect_runtime_error(lambda: read_stdhep(factory, path), "Files containing n-tuples are not supported")


def test_stdhep_too_many_particles(load_delphes, tmp_path):
    _, factory, _ = load_delphes({})
    blob = xdr_u32(101) + xdr_u32(0) + xdr_string("2.00") + xdr_u32(1) + xdr_u32(0xFFFFFFFF)
    path = write_stdhep(tmp_path, "big_event.stdhep", blob)
    expect_runtime_error(lambda: read_stdhep(factory, path), "too many particles in event")


def pileup_merger_config(pileup_file):
    return make_config(
        "PileUpMerger",
        PileUpFile=str(pileup_file),
        MeanPileUp=1,
        PileUpDistribution=2,
        InputArray="Delphes/inputParticles",
    )


def test_pileup_merger_missing_file(load_delphes, tmp_path):
    module, _, _ = load_delphes(pileup_merger_config(tmp_path / "missing.pileup"))
    module.ExportArray("inputParticles")
    expect_runtime_error(module.Init, "can't open pile-up file")


def test_pileup_reader_too_many_events(load_delphes, tmp_path):
    path = tmp_path / "huge_entries.pileup"
    path.write_bytes(xdr_u64(0xFFFFFFFFFFFFFFFF))
    module, _, _ = load_delphes(pileup_merger_config(path))
    module.ExportArray("inputParticles")
    expect_runtime_error(module.Init, "too many events in pile-up file")


def test_pileup_reader_too_many_particles(load_delphes, tmp_path):
    blob = xdr_u32(0xFFFFFFFF) + xdr_u64(0) + xdr_u64(1)
    path = tmp_path / "huge_entry.pileup"
    path.write_bytes(blob)
    module, _, _ = load_delphes(pileup_merger_config(path))
    module.ExportArray("inputParticles")
    module.Init()
    expect_runtime_error(module.Process, "too many particles in pile-up event")
