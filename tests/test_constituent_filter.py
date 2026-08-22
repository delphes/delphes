from conftest import build_config, make_candidate, make_jet


def make_module_config(**extra):
    return build_config(
        "ConstituentFilter",
        {
            "JetPTMin": 0.0,
            "JetInputArray": "Delphes/inputJets",
            "ConstituentInputArray": ["Delphes/inputConstituents", "outputConstituents"],
        },
        **extra,
    )


def run_module_test(run_generic, config, jet_pt=50.0, n_constituents=3):
    def setup(module, factory):
        jet_array = module.ExportArray("inputJets")
        jet = make_jet(factory, jet_pt, 0.5)
        constituent_array = module.ExportArray("inputConstituents")
        for i in range(n_constituents):
            c = make_candidate(factory, 10.0, 0.5, 0.1 * i, pid=211, charge=1)
            constituent_array.Add(c)
            jet.AddCandidate(c)
        jet_array.Add(jet)

    return run_generic(config, setup=setup, outputs=("TestModule/outputConstituents",))


def test_filters_constituents(run_generic):
    output = run_module_test(run_generic, make_module_config(), n_constituents=3)
    assert output.GetEntries() == 3


def test_jet_pt_min_rejects(run_generic):
    output = run_module_test(run_generic, make_module_config(JetPTMin=100.0), jet_pt=50.0, n_constituents=3)
    assert output.GetEntries() == 0


def test_preserves_candidate_properties(run_generic):
    output = run_module_test(run_generic, make_module_config(), n_constituents=1)
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 10.0
