from .conftest import assert_deterministic, candidate_snapshots, make_config


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


def test_unstable_long_trajectory_decays(run_module):
    output = run_module(make_module_config(), [{"pid": 130, "charge": 1, "P": 50.0, "L": 1e12}])
    assert output.GetEntries() == 0


def test_zero_trajectory_never_decays(run_module):
    output = run_module(make_module_config(), [{"pid": 130, "charge": 1, "P": 50.0, "L": 0.0}])
    assert output.GetEntries() == 1


def test_empty_input(run_module):
    output = run_module(make_module_config(), [])
    assert output.GetEntries() == 0


def test_deterministic_with_fixed_seed(run_module):
    particles = [
        {"pid": 130, "charge": 1, "P": 50.0, "L": 1e4},
        {"pid": 130, "charge": 1, "P": 50.0, "L": 3e4},
        {"pid": 130, "charge": 1, "P": 50.0, "L": 1e5},
        {"pid": 130, "charge": 1, "P": 50.0, "L": 3e5},
        {"pid": 130, "charge": 1, "P": 50.0, "L": 1e6},
    ]

    assert_deterministic(
        lambda: run_module(make_module_config(), particles),
        extract=lambda out: candidate_snapshots(out, ("PID", "P")),
    )


def test_zero_pt_unstable_decays_immediately(run_module):
    output = run_module(make_module_config(), [{"pid": 130, "charge": 1, "P": 1e-6, "L": 10.0}])
    assert output.GetEntries() == 0


def test_high_pt_unstable_never_decays(run_module):
    output = run_module(make_module_config(), [{"pid": 130, "charge": 1, "P": 1e10, "L": 10.0}])
    assert output.GetEntries() == 1


def test_extreme_eta_stable_kept(run_module):
    output = run_module(make_module_config(), [{"pid": 22, "charge": 1, "pt": 1e4, "eta": 5.0}])
    assert output.GetEntries() == 1
