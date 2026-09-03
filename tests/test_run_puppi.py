import pytest
from .conftest import assert_deterministic, candidate_snapshots, make_candidate, make_config

PUPPI_PARAMS = {
    "EtaMinBin": [-5.0],
    "EtaMaxBin": [5.0],
    "PtMinBin": [0.0],
    "ConeSizeBin": [0.4],
    "RMSPtMinBin": [0.0],
    "RMSScaleFactorBin": [1.0],
    "NeutralMinEBin": [0.0],
    "NeutralPtSlope": [0.0],
    "UseCharged": [1],
    "ApplyLowPUCorr": [0],
    "MetricId": [0],
    "CombId": [0],
}


def puppi_config(**extra):
    params = dict(PUPPI_PARAMS)
    params.update(extra)
    return make_config(
        "RunPUPPI",
        TrackInputArray="Delphes/inputTracks",
        NeutralInputArray="Delphes/inputNeutrals",
        PVInputArray="Delphes/inputVertices",
        OutputArray="puppiParticles",
        OutputArrayTracks="puppiTracks",
        OutputArrayNeutrals="puppiNeutrals",
        **params,
    )


def add_track(factory, array, pt, is_reco_pu, pid=211):
    c = make_candidate(factory, pt, 0.0, pid=pid, charge=1)
    c.IsRecoPU = is_reco_pu
    constituent = make_candidate(factory, 0.0, 0.0)
    constituent.Position.SetXYZT(0.0, 0.0, 0.0, 0.0)
    c.AddCandidate(constituent)
    array.Add(c)


def add_neutral(factory, array, pt, phi=0.1, pid=22):
    c = make_candidate(factory, pt, 0.0, phi=phi, pid=pid, charge=0)
    constituent = make_candidate(factory, 0.0, 0.0)
    constituent.Position.SetXYZT(0.0, 0.0, 0.0, 0.0)
    c.AddCandidate(constituent)
    array.Add(c)


def run_puppi(load_delphes, config):
    module, factory, _ = load_delphes(config)
    tracks = module.ExportArray("inputTracks")
    add_track(factory, tracks, 30.0, is_reco_pu=0)
    add_track(factory, tracks, 20.0, is_reco_pu=1)
    neutrals = module.ExportArray("inputNeutrals")
    add_neutral(factory, neutrals, 40.0)
    vertices = module.ExportArray("inputVertices")
    vertex = make_candidate(factory, 0.0, 0.0)
    vertex.Position.SetXYZT(0.0, 0.0, 0.0, 0.0)
    vertices.Add(vertex)
    module.Init()
    module.Process()
    return {
        name: module.ImportArray(f"TestModule/{name}") for name in ("puppiParticles", "puppiTracks", "puppiNeutrals")
    }


def test_output_arrays(load_delphes):
    outputs = run_puppi(load_delphes, puppi_config())

    assert outputs["puppiParticles"].GetEntries() == 2
    assert outputs["puppiTracks"].GetEntries() == 1
    assert outputs["puppiNeutrals"].GetEntries() == 1
    assert outputs["puppiTracks"].At(0).Momentum.Pt() == pytest.approx(30.0, rel=1e-6)
    assert outputs["puppiNeutrals"].At(0).Momentum.Pt() == pytest.approx(40.0, rel=1e-6)


def test_pileup_track_is_removed(load_delphes):
    outputs = run_puppi(load_delphes, puppi_config())
    pts = sorted(
        round(outputs["puppiParticles"].At(i).Momentum.Pt(), 6) for i in range(outputs["puppiParticles"].GetEntries())
    )
    assert 30.0 in pts
    assert 40.0 in pts
    assert 20.0 not in pts


def test_min_puppi_weight_threshold(load_delphes):
    outputs = run_puppi(load_delphes, puppi_config(MinPuppiWeight=1.0e9))
    assert outputs["puppiParticles"].GetEntries() == 0
    assert outputs["puppiTracks"].GetEntries() == 0
    assert outputs["puppiNeutrals"].GetEntries() == 0


def test_deterministic_with_fixed_seed(load_delphes):
    config = puppi_config()
    assert_deterministic(
        lambda: run_puppi(load_delphes, config)["puppiParticles"],
        extract=lambda out: candidate_snapshots(out, ("Momentum", "IsPU")),
        abs_tol=1e-9,
    )


def test_use_no_lep_flag(load_delphes):
    outputs = run_puppi(load_delphes, puppi_config(UseNoLep=False))
    assert outputs["puppiParticles"].GetEntries() == 2
