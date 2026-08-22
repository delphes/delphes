from conftest import build_config


def make_module_config(**extra):
    return build_config(
        "JetFlavorAssociation",
        {
            "PartonInputArray": "Delphes/inputPartons",
            "ParticleInputArray": "Delphes/inputParticles",
            "JetInputArray": "Delphes/inputJets",
            "DeltaR": 0.5,
            "PartonPTMin": 0.0,
            "PartonEtaMax": 10.0,
        },
        **extra,
    )


def test_sets_flavor(run_tagging):
    jets = run_tagging(
        make_module_config(),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 0.5, "pid": 5, "status": 3}],
        particles=[],
    )
    assert jets.GetEntries() == 1
    assert jets.At(0).Flavor == 5


def test_c_quark_flavor(run_tagging):
    jets = run_tagging(
        make_module_config(),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 0.5, "pid": 4, "status": 3}],
        particles=[],
    )
    assert jets.At(0).Flavor == 4


def test_no_parton_in_cone(run_tagging):
    jets = run_tagging(
        make_module_config(),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 3.0, "pid": 5, "status": 3}],
        particles=[],
    )
    assert jets.At(0).Flavor == 0


def test_light_quark_flavor(run_tagging):
    jets = run_tagging(
        make_module_config(),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 0.5, "pid": 1, "status": 3}],
        particles=[],
    )
    assert jets.At(0).Flavor == 1
