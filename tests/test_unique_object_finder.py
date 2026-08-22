from conftest import build_config, make_candidate


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
