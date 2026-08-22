import pytest
from conftest import make_config


def test_smearing_zero_resolution(run_module):
    config = make_config("MomentumSmearing", ResolutionFormula=0.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 1
    smeared = output.At(0)
    assert smeared.Momentum.Pt() == pytest.approx(50.0, rel=1e-6)


def test_smearing_preserves_mass(run_module):
    config = make_config("MomentumSmearing", ResolutionFormula=0.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    smeared = output.At(0)
    assert smeared.Momentum.M() == pytest.approx(0.0, abs=1e-6)


def test_smearing_sets_track_resolution(run_module):
    config = make_config("MomentumSmearing", ResolutionFormula=0.1)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    smeared = output.At(0)
    assert smeared.TrackResolution == pytest.approx(0.1, rel=1e-6)


def test_smearing_modifies_pt(run_module):
    config = make_config("MomentumSmearing", ResolutionFormula=0.1)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    smeared = output.At(0)
    assert smeared.Momentum.Pt() != 50.0
