from .conftest import build_config


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


def test_empty_input(run_tagging):
    jets = run_tagging(make_module_config(), jets=[], tracks=[])
    assert jets.GetEntries() == 0


def test_no_b_tag_tracks_far_from_jet(run_tagging):
    config = make_module_config(DeltaR=0.3, SigMin=0.1, Ntracks=1, TrackPtMin=0.0)
    jets = run_tagging(
        config,
        jets=[{"pt": 50.0, "eta": 0.5}],
        tracks=[{"pt": 10.0, "eta": 3.0, "charge": 1, "pid": 211}],
    )
    assert jets.GetEntries() == 1
    assert jets.At(0).BTag & 1 == 0


def tag_track(d0, error_d0):
    return {"pt": 10.0, "eta": 0.5, "charge": 1, "pid": 211, "D0": d0, "ErrorD0": error_d0, "Xd": 1.0, "Yd": 0.0}


def test_b_tag_with_high_significance_track(run_tagging):
    config = make_module_config(Ntracks=1, SigMin=3.0, TrackIPMax=2.0)
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[tag_track(1.9, 0.1)])
    assert jets.GetEntries() == 1
    assert jets.At(0).BTag & 1 == 1


def test_no_b_tag_low_significance_track(run_tagging):
    config = make_module_config(Ntracks=1, SigMin=3.0)
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[tag_track(0.01, 0.1)])
    assert jets.GetEntries() == 1
    assert jets.At(0).BTag & 1 == 0


def test_bit_number_1(run_tagging):
    config = make_module_config(BitNumber=1, Ntracks=1, SigMin=3.0)
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[tag_track(1.9, 0.1)])
    assert jets.At(0).BTag & 1 == 0
    assert jets.At(0).BTag & 2 == 2


def test_use_3d_significance(run_tagging):
    config = make_module_config(Ntracks=1, SigMin=3.0, Use3D=True)
    track = tag_track(1.9, 0.1)
    track["Zd"] = 1.0
    track["DZ"] = 1.0
    track["ErrorDZ"] = 0.1
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[track])
    assert jets.At(0).BTag & 1 == 1


def test_sip_at_sigmin_not_counted(run_tagging):
    config = make_module_config(Ntracks=1, SigMin=2.0)
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[tag_track(1.0, 0.5)])
    assert jets.GetEntries() == 1
    assert jets.At(0).BTag & 1 == 0


def test_sip_just_above_sigmin_counted(run_tagging):
    config = make_module_config(Ntracks=1, SigMin=1.99)
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[tag_track(1.0, 0.5)])
    assert jets.At(0).BTag & 1 == 1


def test_d0_at_ipmax_counted(run_tagging):
    config = make_module_config(Ntracks=1, SigMin=3.0, TrackIPMax=2.0)
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[tag_track(2.0, 0.1)])
    assert jets.At(0).BTag & 1 == 1


def test_d0_above_ipmax_not_counted(run_tagging):
    config = make_module_config(Ntracks=1, SigMin=3.0, TrackIPMax=2.0)
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[tag_track(2.1, 0.1)])
    assert jets.At(0).BTag & 1 == 0


def test_pt_at_track_pt_min_counted(run_tagging):
    track = tag_track(1.9, 0.1)
    track["pt"] = 1.0
    config = make_module_config(Ntracks=1, SigMin=3.0, TrackPtMin=1.0)
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[track])
    assert jets.At(0).BTag & 1 == 1


def test_pt_below_track_pt_min_not_counted(run_tagging):
    track = tag_track(1.9, 0.1)
    track["pt"] = 0.9
    config = make_module_config(Ntracks=1, SigMin=3.0, TrackPtMin=1.0)
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=[track])
    assert jets.At(0).BTag & 1 == 0


def test_count_at_ntracks_inclusive(run_tagging):
    config = make_module_config(Ntracks=3, SigMin=3.0)
    three = [tag_track(1.9, 0.1), tag_track(1.8, 0.2), tag_track(1.7, 0.3)]
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=three)
    assert jets.At(0).BTag & 1 == 1
    jets = run_tagging(config, jets=[{"pt": 50.0, "eta": 0.5}], tracks=three[:2])
    assert jets.At(0).BTag & 1 == 0
