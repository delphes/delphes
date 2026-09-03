from .conftest import build_config, make_candidate, make_jet


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


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputJets")
        module.ExportArray("inputConstituents")

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputConstituents",))
    assert output.GetEntries() == 0


def test_jet_pt_at_boundary_dropped(run_generic):
    output = run_module_test(run_generic, make_module_config(JetPTMin=50.0), jet_pt=50.0, n_constituents=3)
    assert output.GetEntries() == 0


def test_jet_pt_just_above_boundary_kept(run_generic):
    output = run_module_test(run_generic, make_module_config(JetPTMin=49.99), jet_pt=50.0, n_constituents=3)
    assert output.GetEntries() == 3


def run_union_test(run_generic, config, n_jets=2, n_strays=1):
    def setup(module, factory):
        jet_array = module.ExportArray("inputJets")
        constituent_array = module.ExportArray("inputConstituents")
        kept = []
        for j in range(n_jets):
            jet = make_jet(factory, 50.0 + j, 0.5)
            for i in range(3):
                c = make_candidate(factory, 10.0, 0.5, 0.1 * i, pid=211, charge=1)
                constituent_array.Add(c)
                jet.AddCandidate(c)
                kept.append(c)
            jet_array.Add(jet)
        for i in range(n_strays):
            stray = make_candidate(factory, 10.0, 1.0, pid=211, charge=1)
            constituent_array.Add(stray)
            kept.append(stray)

    return run_generic(config, setup=setup, outputs=("TestModule/outputConstituents",))


def test_multiple_jets_union_of_constituents(run_generic):
    output = run_union_test(run_generic, make_module_config(), n_jets=2, n_strays=1)
    assert output.GetEntries() == 6


def test_stray_constituent_flag_left_clear(load_delphes):
    config = make_module_config()
    module, factory, _ = load_delphes(config)
    jet_array = module.ExportArray("inputJets")
    array = module.ExportArray("inputConstituents")
    jet = make_jet(factory, 50.0, 0.5)
    member = make_candidate(factory, 10.0, 0.5, pid=211, charge=1)
    stray = make_candidate(factory, 10.0, 1.0, pid=211, charge=1)
    array.Add(member)
    array.Add(stray)
    jet.AddCandidate(member)
    jet_array.Add(jet)
    module.Init()
    module.Process()
    out = module.ImportArray("TestModule/outputConstituents")
    assert out.GetEntries() == 1
    assert member.IsConstituent == 1
    assert stray.IsConstituent == 0
