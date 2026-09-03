from .conftest import make_config


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


def test_empty_input(run_module):
    output = run_module(make_module_config(), [])
    assert output.GetEntries() == 0


def test_extreme_kinematics_preserved(run_module):
    candidates = [
        {"pt": 1e-6, "eta": 5.0, "pid": 211, "charge": 1, "IsRecoPU": 0},
        {"pt": 1e4, "eta": -5.0, "pid": 13, "charge": -1, "IsRecoPU": 1},
        {"pt": 1e-6, "eta": -5.0, "pid": 11, "charge": -1, "IsRecoPU": 0},
    ]
    output = run_module(make_module_config(), candidates)
    assert output.GetEntries() == 2
    assert output.At(0).Momentum.Pt() == 1e-6
    assert output.At(0).PID == 211
    assert output.At(1).Momentum.Pt() == 1e-6
    assert output.At(1).PID == 11
