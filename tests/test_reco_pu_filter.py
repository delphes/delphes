from conftest import make_config


def make_module_config(**extra):
    return make_config("RecoPuFilter", **extra)


def test_keeps_non_pileup(run_module):
    output = run_module(make_module_config(), [{"pid": 211, "charge": 1, "IsRecoPU": 0}])
    assert output.GetEntries() == 1


def test_rejects_pileup(run_module):
    output = run_module(make_module_config(), [{"pid": 211, "charge": 1, "IsRecoPU": 1}])
    assert output.GetEntries() == 0


def test_mixed_pu_non_pu(run_module):
    candidates = [{"pid": 211, "charge": 1, "IsRecoPU": is_pu} for is_pu in [0, 1, 0, 1, 0]]
    output = run_module(make_module_config(), candidates)
    assert output.GetEntries() == 3
