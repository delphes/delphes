from conftest import build_config, make_candidate, make_jet

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
