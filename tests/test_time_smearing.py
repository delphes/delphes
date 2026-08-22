import pytest
from conftest import make_config


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
