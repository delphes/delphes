from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "Merger",
        {
            "InputArray": ["Delphes/inputA", "Delphes/inputB"],
            "OutputArray": "outputCandidates",
            "MomentumOutputArray": "momentum",
            "EnergyOutputArray": "energy",
        },
        **extra,
    )


def run_merger_test(run_generic, config, arrays_data):
    def setup(module, factory):
        for array_name, candidates in arrays_data:
            arr = module.ExportArray(array_name)
            for pt, eta in candidates:
                c = make_candidate(factory, pt, eta)
                arr.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/outputCandidates",))


def test_merge_two_arrays(run_generic):
    arrays = [("inputA", [(50.0, 0.5)]), ("inputB", [(30.0, 1.0)])]
    output = run_merger_test(run_generic, make_module_config(), arrays)
    assert output.GetEntries() == 2


def test_merge_single_array(run_generic):
    arrays = [("inputA", [(50.0, 0.5), (30.0, 1.0)])]
    output = run_merger_test(run_generic, make_module_config(InputArray=["Delphes/inputA"]), arrays)
    assert output.GetEntries() == 2


def test_empty_input(run_generic):
    arrays = [("inputA", []), ("inputB", [])]
    output = run_merger_test(run_generic, make_module_config(), arrays)
    assert output.GetEntries() == 0
