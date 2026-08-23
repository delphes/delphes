from conftest import add_input_jets, assert_deterministic, build_config, candidate_snapshots


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


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputJets")

    results = run_generic(
        make_module_config(),
        setup=setup,
        outputs={
            "jets": "TestModule/jets",
            "electrons": "TestModule/fakeElectrons",
            "muons": "TestModule/fakeMuons",
            "photons": "TestModule/fakePhotons",
        },
    )
    assert results["jets"].GetEntries() == 0
    assert results["electrons"].GetEntries() == 0
    assert results["muons"].GetEntries() == 0
    assert results["photons"].GetEntries() == 0


def test_deterministic_with_fixed_seed(run_generic):
    config = make_module_config(EfficiencyFormula=[11, 0.4, 13, 0.3, 22, 0.2])

    def setup(module, factory):
        add_input_jets(module, factory, [{"pt": 50.0, "eta": 0.5}, {"pt": 60.0, "eta": 0.6}, {"pt": 70.0, "eta": 0.7}])

    outputs = {
        "jets": "TestModule/jets",
        "electrons": "TestModule/fakeElectrons",
        "muons": "TestModule/fakeMuons",
        "photons": "TestModule/fakePhotons",
    }

    def snap(results):
        out = []
        for key in ("jets", "electrons", "muons", "photons"):
            out.extend((key, *row) for row in candidate_snapshots(results[key], ("Momentum", "Charge")))
        return tuple(out)

    assert_deterministic(lambda: run_generic(config, setup=setup, outputs=outputs), extract=snap)
