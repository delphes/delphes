import pytest
from conftest import C_LIGHT, assert_deterministic, candidate_snapshots, make_config, mean, run_repeated, stddev


def test_zero_resolution(run_module):
    config = make_config("TimeSmearing", TimeResolution=0.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 1
    smeared = output.At(0)
    assert smeared.Position.T() == pytest.approx(0.0, abs=1e-6)


def test_smears_time(run_module):
    config = make_config("TimeSmearing", TimeResolution=1.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    smeared = output.At(0)
    assert smeared.Position.T() != 0.0


def test_sets_error_t(run_module):
    config = make_config("TimeSmearing", TimeResolution=1.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    smeared = output.At(0)
    assert smeared.ErrorT > 0


def test_preserves_momentum(run_module):
    config = make_config("TimeSmearing", TimeResolution=1.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    smeared = output.At(0)
    assert smeared.Momentum.Pt() == pytest.approx(50.0, rel=1e-6)


def test_preserves_pid(run_module):
    config = make_config("TimeSmearing", TimeResolution=1.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5, "pid": 211, "charge": 1}])
    assert output.At(0).PID == 211


def test_statistical_mean_and_rms(load_delphes):
    config = make_config("TimeSmearing", TimeResolution=1.0)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = factory.NewCandidate()
        c.Momentum.SetPtEtaPhiE(50.0, 0.5, 0.0, 50.0 * 1.5)
        c.Position.SetXYZT(0.0, 0.0, 0.0, 0.0)
        input_array.Add(c)

    results = run_repeated(load_delphes, config, 200, setup)
    times = [results[i].At(0).Position.T() for i in range(200)]
    assert mean(times) == pytest.approx(0.0, abs=5.0e10)
    assert stddev(times) == pytest.approx(1.0e3 * C_LIGHT, rel=0.1)


def test_deterministic_with_fixed_seed(run_module):
    config = make_config("TimeSmearing", TimeResolution=1.0)
    candidates = [{"pt": 50.0, "eta": 0.5}]
    assert_deterministic(
        lambda: run_module(config, candidates),
        extract=lambda out: candidate_snapshots(out, ("Position", "ErrorT")),
        abs_tol=1e-12,
    )


def test_empty_input(run_module):
    output = run_module(make_config("TimeSmearing", TimeResolution=1.0), [])
    assert output.GetEntries() == 0
