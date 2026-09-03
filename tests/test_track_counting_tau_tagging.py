from .conftest import build_config


def make_module_config(**extra):
    return build_config(
        "TrackCountingTauTagging",
        {
            "ParticleInputArray": "Delphes/inputParticles",
            "PartonInputArray": "Delphes/inputPartons",
            "TrackInputArray": "Delphes/inputTracks",
            "JetInputArray": "Delphes/inputJets",
            "BitNumber": 0,
            "DeltaR": 0.5,
            "DeltaRTrack": 0.3,
            "TrackPTMin": 0.0,
            "TauPTMin": 0.0,
            "TauEtaMax": 10.0,
            "EfficiencyFormula": [[15, 1.0]],
        },
        **extra,
    )


def test_tau_tagging(run_tagging):
    jets = run_tagging(
        make_module_config(),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}],
        particles=[
            {"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1},
            {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1},
        ],
        tracks=[{"pt": 10.0, "eta": 0.5, "charge": 1, "pid": 211}, {"pt": 10.0, "eta": 0.5, "charge": 1, "pid": 211}],
    )
    assert jets.GetEntries() == 1


def test_empty_input(run_tagging):
    jets = run_tagging(make_module_config(), jets=[], partons=[], particles=[], tracks=[])
    assert jets.GetEntries() == 0


def make_tau(parton_eta=0.5):
    return {"pt": 50.0, "eta": parton_eta, "pid": 15, "status": 3, "d1": 0, "d2": 1}


def test_identifier_sign_flip(run_tagging):
    pions = [{"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1}]
    track = {"pt": 10.0, "eta": 0.5, "charge": 1, "pid": 211}
    config = make_module_config(EfficiencyFormula=[1, 0.0, -1, 1.0, 0, 0.0])

    jets = run_tagging(
        config,
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[make_tau(0.5)],
        particles=pions,
        tracks=[track],
    )
    assert jets.At(0).TauTag & 1 == 0

    far_pions = [{"pt": 20.0, "eta": 2.9, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 3.1, "pid": -211, "status": 1}]
    jets = run_tagging(
        config,
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[make_tau(3.0)],
        particles=far_pions,
        tracks=[track],
    )
    assert jets.At(0).TauTag & 1 == 1


def test_identifier_clamped_at_plus_two(run_tagging):
    pions = [{"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1}]
    tracks = [
        {"pt": 10.0, "eta": 0.5, "phi": 0.0, "charge": 1, "pid": 211},
        {"pt": 10.0, "eta": 0.5, "phi": 0.1, "charge": 1, "pid": 211},
        {"pt": 10.0, "eta": 0.5, "phi": 0.2, "charge": 1, "pid": 211},
    ]
    config = make_module_config(EfficiencyFormula=[2, 1.0, 0, 0.0])
    jets = run_tagging(
        config,
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[make_tau(0.5)],
        particles=pions,
        tracks=tracks,
    )
    assert jets.At(0).TauTag & 1 == 1


def test_jet_charge_sum_of_matched_tracks(run_tagging):
    pions = [{"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1}]
    tracks = [
        {"pt": 10.0, "eta": 0.5, "phi": 0.0, "charge": 1, "pid": 211},
        {"pt": 10.0, "eta": 0.5, "phi": 0.1, "charge": -1, "pid": 211},
    ]
    jets = run_tagging(
        make_module_config(),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[make_tau(0.5)],
        particles=pions,
        tracks=tracks,
    )
    assert jets.At(0).Charge == 0


def test_jet_charge_excludes_low_pt_tracks(run_tagging):
    pions = [{"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1}]
    tracks = [
        {"pt": 10.0, "eta": 0.5, "phi": 0.0, "charge": 1, "pid": 211},
        {"pt": 4.0, "eta": 0.5, "phi": 0.1, "charge": 1, "pid": 211},
    ]
    jets = run_tagging(
        make_module_config(TrackPTMin=5.0),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[make_tau(0.5)],
        particles=pions,
        tracks=tracks,
    )
    assert jets.At(0).Charge == 1
