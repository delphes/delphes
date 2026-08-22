from conftest import build_config


def make_module_config(**extra):
    return build_config("IdentificationMap", {"EfficiencyFormula": [13, 11, 1.0]}, **extra)


def test_remapping(run_module):
    output = run_module(make_module_config(), [{"pid": 13, "charge": 1}])
    assert output.GetEntries() == 1
    assert output.At(0).PID == 11


def test_no_match_passes_through(run_module):
    output = run_module(make_module_config(), [{"pid": 211, "charge": 1}])
    assert output.GetEntries() == 1
    assert output.At(0).PID == 211


def test_zero_efficiency_rejects(run_module):
    output = run_module(make_module_config(EfficiencyFormula=[13, 11, 0.0]), [{"pid": 13, "charge": 1}])
    assert output.GetEntries() == 0
