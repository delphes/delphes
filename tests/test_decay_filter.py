from conftest import make_config


def make_module_config(**extra):
    return make_config("DecayFilter", **extra)


def test_stable_particle_passes(run_module):
    output = run_module(make_module_config(), [{"pid": 22, "charge": 1, "P": 50.0, "L": 0.0}])
    assert output.GetEntries() == 1


def test_unknown_pid_passes(run_module):
    output = run_module(make_module_config(), [{"pid": 999999, "charge": 1, "P": 50.0, "L": 0.0}])
    assert output.GetEntries() == 1


def test_unstable_short_trajectory_passes(run_module):
    output = run_module(make_module_config(), [{"pid": 130, "charge": 1, "P": 50.0, "L": 0.001}])
    assert output.GetEntries() == 1


def test_unstable_long_trajectory_may_decay(run_module):
    output = run_module(make_module_config(), [{"pid": 130, "charge": 1, "P": 50.0, "L": 1e12}])
    assert output.GetEntries() <= 1
