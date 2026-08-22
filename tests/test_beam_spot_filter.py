from conftest import make_config

PU = {"pid": 211, "charge": 1, "IsPU": 1}
NON_PU = {"pid": 211, "charge": 1, "IsPU": 0}


def make_module_config(**extra):
    return make_config("BeamSpotFilter", **extra)


def test_passes_non_pu(run_module):
    output = run_module(make_module_config(), [NON_PU])
    assert output.GetEntries() == 1


def test_passes_until_first_non_pu(run_module):
    output = run_module(make_module_config(), [PU, PU, NON_PU])
    assert output.GetEntries() == 3


def test_all_pu(run_module):
    output = run_module(make_module_config(), [PU] * 3)
    assert output.GetEntries() == 3


def test_stops_after_first_non_pu(run_module):
    output = run_module(make_module_config(), [NON_PU] * 3)
    assert output.GetEntries() == 1


def test_preserves_particle_properties(run_module):
    output = run_module(make_module_config(), [{"pt": 42.0, "eta": 1.5, "pid": 13, "charge": -1, "IsPU": 0}])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 42.0
    assert output.At(0).PID == 13
