from .conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "UniqueObjectFinder",
        {
            "UseUniqueID": True,
            "InputArray": ["Delphes/inputElectrons", "outputElectrons"],
        },
        **extra,
    )


def test_unique_objects_pass_through(run_generic):
    config = make_module_config(
        InputArray=["Delphes/inputElectrons", "outputElectrons", "Delphes/inputPhotons", "outputPhotons"]
    )

    def setup(module, factory):
        arr1 = module.ExportArray("inputElectrons")
        c1 = make_candidate(factory, 50.0, 0.5, pid=11, charge=-1)
        arr1.Add(c1)

        arr2 = module.ExportArray("inputPhotons")
        c2 = make_candidate(factory, 30.0, 1.5, pid=22, charge=0)
        arr2.Add(c2)

    electrons, photons = run_generic(
        config, setup=setup, outputs=("TestModule/outputElectrons", "TestModule/outputPhotons")
    )
    assert electrons.GetEntries() == 1
    assert photons.GetEntries() == 1


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputElectrons")

    electrons = run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputElectrons",))
    assert electrons.GetEntries() == 0


def test_preserves_particle(run_generic):
    def setup(module, factory):
        arr = module.ExportArray("inputElectrons")
        c = make_candidate(factory, 42.0, 1.0, pid=11, charge=-1)
        arr.Add(c)

    electrons = run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputElectrons",))
    assert electrons.GetEntries() == 1
    assert electrons.At(0).Momentum.Pt() == 42.0


def run_two_jet_test(run_generic, use_unique_id, share_constituent):
    config = build_config(
        "UniqueObjectFinder",
        {
            "UseUniqueID": use_unique_id,
            "InputArray": [
                "Delphes/inputJets1",
                "outJets1",
                "Delphes/inputJets2",
                "outJets2",
            ],
        },
    )

    def setup(module, factory):
        arr1 = module.ExportArray("inputJets1")
        arr2 = module.ExportArray("inputJets2")
        const1 = make_candidate(factory, 5.0, 0.5, pid=211, charge=1)
        jet1 = make_candidate(factory, 50.0, 0.5, pid=0)
        jet1.AddCandidate(const1)
        arr1.Add(jet1)

        if share_constituent:
            const2 = const1
        else:
            const2 = make_candidate(factory, 5.0, 0.5, pid=211, charge=1)
        jet2 = make_candidate(factory, 40.0, 0.5, pid=0)
        jet2.AddCandidate(const2)
        arr2.Add(jet2)

    out1, out2 = run_generic(config, setup=setup, outputs=("TestModule/outJets1", "TestModule/outJets2"))
    return out1, out2


def test_default_dedup_drops_shared_constituent(run_generic):
    out1, out2 = run_two_jet_test(run_generic, use_unique_id=False, share_constituent=True)
    assert out1.GetEntries() == 1
    assert out2.GetEntries() == 0


def test_default_dedup_keeps_distinct_jets(run_generic):
    out1, out2 = run_two_jet_test(run_generic, use_unique_id=False, share_constituent=False)
    assert out1.GetEntries() == 1
    assert out2.GetEntries() == 1
    assert out2.At(0).Momentum.Pt() == 40.0


def test_use_unique_id_ignores_shared_constituent(run_generic):
    out1, out2 = run_two_jet_test(run_generic, use_unique_id=True, share_constituent=True)
    assert out1.GetEntries() == 1
    assert out2.GetEntries() == 1
    assert out2.At(0).Momentum.Pt() == 40.0


def test_same_candidate_in_both_arrays_dropped_from_second(run_generic):
    config = build_config(
        "UniqueObjectFinder",
        {
            "UseUniqueID": False,
            "InputArray": [
                "Delphes/inputElectrons",
                "outputElectrons",
                "Delphes/inputPhotons",
                "outputPhotons",
            ],
        },
    )

    def setup(module, factory):
        arr1 = module.ExportArray("inputElectrons")
        arr2 = module.ExportArray("inputPhotons")
        c = make_candidate(factory, 50.0, 0.5, pid=11, charge=-1)
        arr1.Add(c)
        arr2.Add(c)

    electrons, photons = run_generic(
        config, setup=setup, outputs=("TestModule/outputElectrons", "TestModule/outputPhotons")
    )
    assert electrons.GetEntries() == 1
    assert electrons.At(0).Momentum.Pt() == 50.0
    assert photons.GetEntries() == 0


def test_extreme_kinematics_distinct_candidates_kept(run_generic):
    config = make_module_config(
        InputArray=["Delphes/inputElectrons", "outputElectrons", "Delphes/inputPhotons", "outputPhotons"]
    )

    def setup(module, factory):
        arr1 = module.ExportArray("inputElectrons")
        arr1.Add(make_candidate(factory, 1e-6, 5.0, pid=11, charge=-1))
        arr1.Add(make_candidate(factory, 1e4, -5.0, pid=11, charge=-1))
        arr2 = module.ExportArray("inputPhotons")
        arr2.Add(make_candidate(factory, 1e-6, 5.0, pid=22, charge=0))

    electrons, photons = run_generic(
        config, setup=setup, outputs=("TestModule/outputElectrons", "TestModule/outputPhotons")
    )
    assert electrons.GetEntries() == 2
    assert electrons.At(0).Momentum.Pt() == 1e-6
    assert electrons.At(1).Momentum.Pt() == 1e4
    assert photons.GetEntries() == 1
    assert photons.At(0).Momentum.Pt() == 1e-6


def test_extreme_kinematics_shared_constituent_dropped(run_generic):
    config = build_config(
        "UniqueObjectFinder",
        {
            "UseUniqueID": False,
            "InputArray": ["Delphes/inputJets1", "outJets1", "Delphes/inputJets2", "outJets2"],
        },
    )

    def setup(module, factory):
        arr1 = module.ExportArray("inputJets1")
        arr2 = module.ExportArray("inputJets2")
        const1 = make_candidate(factory, 5.0, 0.5, pid=211, charge=1)
        jet1 = make_candidate(factory, 1e4, 5.0, pid=0)
        jet1.AddCandidate(const1)
        arr1.Add(jet1)

        jet2 = make_candidate(factory, 1e-6, -5.0, pid=0)
        jet2.AddCandidate(const1)
        arr2.Add(jet2)

    out1, out2 = run_generic(config, setup=setup, outputs=("TestModule/outJets1", "TestModule/outJets2"))
    assert out1.GetEntries() == 1
    assert out2.GetEntries() == 0
