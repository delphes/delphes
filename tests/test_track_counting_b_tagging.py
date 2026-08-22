from conftest import build_config


def make_module_config(**extra):
    return build_config(
        "TrackCountingBTagging",
        {
            "TrackInputArray": "Delphes/inputTracks",
            "JetInputArray": "Delphes/inputJets",
            "BitNumber": 0,
            "TrackPtMin": 1.0,
            "DeltaR": 0.5,
            "TrackIPMax": 2.0,
            "SigMin": 3.0,
            "Ntracks": 3,
            "Use3D": False,
        },
        **extra,
    )


def test_no_b_tag_without_enough_tracks(run_tagging):
    jets = run_tagging(
        make_module_config(),
        jets=[{"pt": 50.0, "eta": 0.5}],
        tracks=[{"pt": 10.0, "eta": 0.5, "charge": 1, "pid": 211}],
    )
    assert jets.GetEntries() == 1
    assert jets.At(0).BTag & 1 == 0


def test_no_b_tag_tracks_far_from_jet(run_tagging):
    config = make_module_config(DeltaR=0.3, SigMin=0.1, Ntracks=1, TrackPtMin=0.0)
    jets = run_tagging(
        config,
        jets=[{"pt": 50.0, "eta": 0.5}],
        tracks=[{"pt": 10.0, "eta": 3.0, "charge": 1, "pid": 211}],
    )
    assert jets.GetEntries() == 1
    assert jets.At(0).BTag & 1 == 0
