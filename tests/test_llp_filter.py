from conftest import build_config


def make_module_config(**extra):
    return build_config("LLPFilter", {"PdgCode": [1000022]}, **extra)


def test_keeps_matching_llp(run_module):
    output = run_module(make_module_config(), [{"pid": 1000022, "charge": 1, "D1": -1, "D2": -1}])
    assert output.GetEntries() == 1


def test_rejects_non_matching(run_module):
    output = run_module(make_module_config(), [{"pid": 13, "charge": 1, "D1": -1, "D2": -1}])
    assert output.GetEntries() == 0


def test_pt_min_filter(run_module):
    output = run_module(
        make_module_config(PTMin=100.0), [{"pid": 1000022, "charge": 1, "D1": -1, "D2": -1, "pt": 50.0}]
    )
    assert output.GetEntries() == 0


def test_daughter_number_filter(run_module):
    output = run_module(make_module_config(DaughterNumber=2), [{"pid": 1000022, "charge": 1, "D1": -1, "D2": -1}])
    assert output.GetEntries() == 0
