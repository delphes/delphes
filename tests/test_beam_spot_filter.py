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


def test_empty_input(run_module):
    output = run_module(make_module_config(), [])
    assert output.GetEntries() == 0


def test_extreme_kinematics_flag_logic_unchanged(run_module):
    output = run_module(
        make_module_config(),
        [
            {"pt": 1e-6, "eta": 5.0, "pid": 211, "charge": 1, "IsPU": 1},
            {"pt": 1e4, "eta": -5.0, "pid": 13, "charge": -1, "IsPU": 0},
        ],
    )
    assert output.GetEntries() == 2
    assert output.At(0).Momentum.Pt() == 1e-6
    assert output.At(0).IsPU == 1
    assert output.At(1).Momentum.Pt() == 1e4
    assert output.At(1).PID == 13
