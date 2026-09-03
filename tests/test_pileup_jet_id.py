import math

import pytest
from .conftest import build_config, make_candidate, make_jet

PILEUP_JET_ID_DEFAULTS = {
    "JetInputArray": "Delphes/inputJets",
    "TrackInputArray": "Delphes/inputTracks",
    "NeutralInputArray": "Delphes/inputNeutrals",
    "OutputArray": "jets",
    "NeutralsInPassingJets": "eflowtowers",
    "JetPTMin": 0.0,
    "UseConstituents": 0,
    "MeanSqDeltaRMaxBarrel": 10.0,
    "BetaMinBarrel": 0.0,
    "MeanSqDeltaRMaxEndcap": 10.0,
    "BetaMinEndcap": 0.0,
    "MeanSqDeltaRMaxForward": 10.0,
    "JetPTMinForNeutrals": 0.0,
    "NeutralPTMin": 0.0,
}


def make_module_config(**extra):
    return build_config("PileUpJetID", PILEUP_JET_ID_DEFAULTS, **extra)


def run_pujetid_test(run_generic, config, jet_pt=50.0, n_tracks=5, track_is_pu_list=None):
    if track_is_pu_list is None:
        track_is_pu_list = [0] * n_tracks

    def setup(module, factory):
        jet_array = module.ExportArray("inputJets")
        jet = make_jet(factory, jet_pt, 0.5)

        track_array = module.ExportArray("inputTracks")
        for i, is_pu in enumerate(track_is_pu_list):
            t = make_candidate(factory, 10.0, 0.5, 0.1 * i, pid=211, charge=1)
            t.IsPU = is_pu
            track_array.Add(t)
            jet.AddCandidate(t)

        jet_array.Add(jet)
        module.ExportArray("inputNeutrals")

    return run_generic(config, setup=setup, outputs=("TestModule/jets",))


def test_jet_passes_puid(run_generic):
    output = run_pujetid_test(run_generic, make_module_config(), n_tracks=5, track_is_pu_list=[0, 0, 0, 0, 0])
    assert output.GetEntries() == 1


def run_pujetid_full(
    run_generic,
    config,
    jet_pt=50.0,
    n_tracks=5,
    track_is_pu_list=None,
    neutral_pt=None,
    use_constituents=False,
    phi_step=0.1,
):
    if track_is_pu_list is None:
        track_is_pu_list = [0] * n_tracks

    def setup(module, factory):
        jet_array = module.ExportArray("inputJets")
        jet = make_jet(factory, jet_pt, 0.5)

        track_array = module.ExportArray("inputTracks")
        for i, is_pu in enumerate(track_is_pu_list):
            t = make_candidate(factory, 10.0, 0.5, phi_step * i, pid=211, charge=1)
            t.IsPU = is_pu
            t.IsRecoPU = is_pu
            track_array.Add(t)
            if use_constituents:
                jet.AddCandidate(t)

        neutral_array = module.ExportArray("inputNeutrals")
        if neutral_pt is not None:
            n = make_candidate(factory, neutral_pt, 0.5, 0.0, pid=22, charge=0)
            neutral_array.Add(n)
            if use_constituents:
                jet.AddCandidate(n)

        jet_array.Add(jet)

    return run_generic(config, setup=setup, outputs={"jets": "TestModule/jets", "eflow": "TestModule/eflowtowers"})


def test_rejects_pileup_dominated_jet(run_generic):
    config = make_module_config(BetaMinBarrel=0.5, MeanSqDeltaRMaxBarrel=10.0)
    results = run_pujetid_full(run_generic, config, n_tracks=5, track_is_pu_list=[1] * 5)
    assert results["jets"].GetEntries() == 1
    assert results["jets"].At(0).Beta == pytest.approx(0.0, abs=1e-6)

    assert results["eflow"].GetEntries() == 0


def test_passes_with_pv_tracks_and_neutral(run_generic):
    config = make_module_config(BetaMinBarrel=0.5, MeanSqDeltaRMaxBarrel=10.0)
    results = run_pujetid_full(run_generic, config, n_tracks=5, track_is_pu_list=[0] * 5, neutral_pt=5.0)
    assert results["jets"].At(0).Beta == pytest.approx(1.0, abs=1e-6)
    assert results["eflow"].GetEntries() == 1


def test_rejects_high_mean_sq_delta_r(run_generic):
    config = make_module_config(BetaMinBarrel=0.5, MeanSqDeltaRMaxBarrel=0.01)
    results = run_pujetid_full(run_generic, config, n_tracks=5, track_is_pu_list=[0] * 5, neutral_pt=5.0)

    assert results["eflow"].GetEntries() == 0


def test_use_constituents_neutral_output(run_generic):
    config = make_module_config(
        UseConstituents=1, BetaMinBarrel=0.5, MeanSqDeltaRMaxBarrel=10.0, JetPTMinForNeutrals=0.0, NeutralPTMin=0.0
    )
    results = run_pujetid_full(
        run_generic, config, n_tracks=5, track_is_pu_list=[0] * 5, neutral_pt=5.0, use_constituents=True
    )
    assert results["jets"].At(0).Beta == pytest.approx(1.0, abs=1e-6)
    assert results["eflow"].GetEntries() == 1


def run_pujetid_boundary_test(run_generic, config, jet_pt=50.0, jet_eta=0.5, tracks=(), neutral_pt=5.0):
    use_constituents = config.get("TestModule", {}).get("UseConstituents")

    def setup(module, factory):
        jet_array = module.ExportArray("inputJets")
        jet = make_jet(factory, jet_pt, jet_eta)
        track_array = module.ExportArray("inputTracks")
        for pt, phi, is_pu in tracks:
            t = make_candidate(factory, pt, jet_eta, phi, pid=211, charge=1)
            t.IsPU = is_pu
            t.IsRecoPU = is_pu
            track_array.Add(t)
            if use_constituents:
                jet.AddCandidate(t)
        neutral_array = module.ExportArray("inputNeutrals")
        if neutral_pt is not None:
            n = make_candidate(factory, neutral_pt, jet_eta, 0.0, pid=22, charge=0)
            neutral_array.Add(n)
            if use_constituents:
                jet.AddCandidate(n)
        jet_array.Add(jet)

    return run_generic(config, setup=setup, outputs={"jets": "TestModule/jets", "eflow": "TestModule/eflowtowers"})


def test_beta_at_min_strict(run_generic):
    tracks = [(10.0, 0.0, 0), (10.0, 0.2, 1)]
    at_cut = run_pujetid_boundary_test(
        run_generic, make_module_config(BetaMinBarrel=0.5, MeanSqDeltaRMaxBarrel=10.0), tracks=tracks
    )
    assert at_cut["jets"].At(0).Beta == pytest.approx(0.5, abs=1e-6)
    assert at_cut["eflow"].GetEntries() == 0
    below = run_pujetid_boundary_test(
        run_generic, make_module_config(BetaMinBarrel=0.49, MeanSqDeltaRMaxBarrel=10.0), tracks=tracks
    )
    assert below["eflow"].GetEntries() == 1


def test_mean_sq_delta_r_at_max_strict(run_generic):
    tracks = [(10.0, 0.0, 0)]
    at_cut = run_pujetid_boundary_test(
        run_generic, make_module_config(BetaMinBarrel=0.0, MeanSqDeltaRMaxBarrel=0.0), tracks=tracks
    )
    assert at_cut["jets"].At(0).MeanSqDeltaR == 0.0
    assert at_cut["eflow"].GetEntries() == 0
    below = run_pujetid_boundary_test(
        run_generic, make_module_config(BetaMinBarrel=0.0, MeanSqDeltaRMaxBarrel=0.01), tracks=tracks
    )
    assert below["eflow"].GetEntries() == 1


def test_eta_region_barrel_endcap_boundary(run_generic):
    tracks = [(10.0, 0.0, 0), (10.0, 0.2, 1)]
    config = make_module_config(
        BetaMinBarrel=0.9, BetaMinEndcap=0.0, MeanSqDeltaRMaxBarrel=10.0, MeanSqDeltaRMaxEndcap=10.0
    )
    barrel = run_pujetid_boundary_test(run_generic, config, jet_eta=1.49, tracks=tracks)
    assert barrel["eflow"].GetEntries() == 0
    endcap = run_pujetid_boundary_test(run_generic, config, jet_eta=1.51, tracks=tracks)
    assert endcap["eflow"].GetEntries() == 1


def test_eta_region_endcap_forward_boundary(run_generic):
    tracks = [(10.0, 0.0, 0), (10.0, 0.2, 1)]
    config = make_module_config(
        BetaMinBarrel=0.0, BetaMinEndcap=0.0, MeanSqDeltaRMaxEndcap=10.0, MeanSqDeltaRMaxForward=0.0
    )
    endcap = run_pujetid_boundary_test(run_generic, config, jet_eta=3.99, tracks=tracks)
    assert endcap["eflow"].GetEntries() == 1
    forward = run_pujetid_boundary_test(run_generic, config, jet_eta=4.01, tracks=tracks)
    assert forward["eflow"].GetEntries() == 0


def test_jet_pt_at_min_for_neutrals(run_generic):
    tracks = [(10.0, 0.0, 0)]
    at_cut = run_pujetid_boundary_test(
        run_generic,
        make_module_config(BetaMinBarrel=0.0, MeanSqDeltaRMaxBarrel=10.0, JetPTMinForNeutrals=50.0),
        jet_pt=50.0,
        tracks=tracks,
    )
    assert at_cut["eflow"].GetEntries() == 0
    above = run_pujetid_boundary_test(
        run_generic,
        make_module_config(BetaMinBarrel=0.0, MeanSqDeltaRMaxBarrel=10.0, JetPTMinForNeutrals=49.9),
        jet_pt=50.0,
        tracks=tracks,
    )
    assert above["eflow"].GetEntries() == 1


def test_fracpt_annulus_bins(run_generic):
    config = make_module_config(UseConstituents=1)
    results = run_pujetid_boundary_test(run_generic, config, tracks=[(10.0, 0.0, 0)], neutral_pt=None)
    jet = results["jets"].At(0)
    assert jet.PTD == pytest.approx(1.0, abs=1e-6)
    for i in range(5):
        assert jet.FracPt[i] == 0.0

    results = run_pujetid_boundary_test(run_generic, config, tracks=[(10.0, 0.05, 0), (10.0, 0.15, 0)], neutral_pt=None)
    jet = results["jets"].At(0)
    assert jet.FracPt[0] == pytest.approx(0.5, abs=1e-6)
    assert jet.FracPt[1] == pytest.approx(0.5, abs=1e-6)
    for i in range(2, 5):
        assert jet.FracPt[i] == 0.0
    assert jet.PTD == pytest.approx(math.sqrt(200.0) / 20.0, rel=1e-5)
