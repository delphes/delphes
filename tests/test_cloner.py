import pytest
from .conftest import make_config


def test_clones_single_candidate(run_module):
    config = make_config("Cloner")
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 50.0


def test_clones_multiple_candidates(run_module):
    config = make_config("Cloner")
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}, {"pt": 30.0, "eta": 1.0}, {"pt": 10.0, "eta": 2.0}])
    assert output.GetEntries() == 3
    assert output.At(0).Momentum.Pt() == 50.0
    assert output.At(1).Momentum.Pt() == 30.0
    assert output.At(2).Momentum.Pt() == 10.0


def test_preserves_eta_phi(run_module):
    config = make_config("Cloner")
    output = run_module(config, [{"pt": 50.0, "eta": 1.2, "phi": 0.7}])
    smeared = output.At(0)
    assert smeared.Momentum.Eta() == pytest.approx(1.2, rel=1e-6)
    assert smeared.Momentum.Phi() == pytest.approx(0.7, rel=1e-6)


def test_empty_input(run_module):
    output = run_module(make_config("Cloner"), [])
    assert output.GetEntries() == 0
