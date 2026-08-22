from conftest import add_input_jets, build_config


def make_module_config(**extra):
    return build_config(
        "JetFakeParticle",
        {
            "InputArray": "Delphes/inputJets",
            "JetOutputArray": "jets",
            "ElectronOutputArray": "fakeElectrons",
            "MuonOutputArray": "fakeMuons",
            "PhotonOutputArray": "fakePhotons",
            "EfficiencyFormula": [11, 0.0, 13, 0.0, 22, 0.0],
        },
        **extra,
    )


def run_fake_test(run_generic, config, jet_pt=50.0, jet_eta=0.5):
    def setup(module, factory):
        add_input_jets(module, factory, [{"pt": jet_pt, "eta": jet_eta}])

    return run_generic(
        config,
        setup=setup,
        outputs={
            "jets": "TestModule/jets",
            "electrons": "TestModule/fakeElectrons",
            "muons": "TestModule/fakeMuons",
            "photons": "TestModule/fakePhotons",
        },
    )


def test_jet_passes_through(run_generic):
    results = run_fake_test(run_generic, make_module_config())
    assert results["jets"].GetEntries() == 1
    assert results["electrons"].GetEntries() == 0
    assert results["muons"].GetEntries() == 0
    assert results["photons"].GetEntries() == 0


def test_fake_electrons_generated(run_generic):
    results = run_fake_test(run_generic, make_module_config(EfficiencyFormula=[11, 1.0]))
    assert results["jets"].GetEntries() == 0
    assert results["electrons"].GetEntries() == 1


def test_fake_muons_generated(run_generic):
    results = run_fake_test(run_generic, make_module_config(EfficiencyFormula=[13, 1.0]))
    assert results["jets"].GetEntries() == 0
    assert results["muons"].GetEntries() == 1


def test_fake_photons_generated(run_generic):
    results = run_fake_test(run_generic, make_module_config(EfficiencyFormula=[22, 1.0]))
    assert results["jets"].GetEntries() == 0
    assert results["photons"].GetEntries() == 1
