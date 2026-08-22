from conftest import build_config


def make_module_config(**extra):
    return build_config(
        "TauTagging",
        {
            "ParticleInputArray": "Delphes/inputParticles",
            "PartonInputArray": "Delphes/inputPartons",
            "JetInputArray": "Delphes/inputJets",
            "BitNumber": 0,
            "DeltaR": 0.5,
            "TauPTMin": 0.0,
            "TauEtaMax": 10.0,
            "EfficiencyFormula": [15, 1.0],
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
    )
    assert jets.GetEntries() == 1
    assert jets.At(0).TauTag & 1 == 1


def test_no_tau_tag_light_jet(run_tagging):
    jets = run_tagging(
        make_module_config(EfficiencyFormula=[15, 1.0, 0, 0.0]),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 0.5, "pid": 1, "status": 3, "d1": -1, "d2": -1}],
        particles=[],
    )
    assert jets.GetEntries() == 1
    assert jets.At(0).TauTag & 1 == 0
