import pytest
from conftest import make_config


def test_zero_resolution(run_module):
    config = make_config("AngularSmearing", EtaResolutionFormula=0.0, PhiResolutionFormula=0.0)
    output = run_module(config, [{"pt": 50.0, "eta": 1.2}])
    smeared = output.At(0)
    assert smeared.Momentum.Eta() == pytest.approx(1.2, rel=1e-6)


def test_smears_eta(run_module):
    config = make_config("AngularSmearing", EtaResolutionFormula=0.05, PhiResolutionFormula=0.0)
    output = run_module(config, [{"pt": 50.0, "eta": 1.2}])
    smeared = output.At(0)
    assert smeared.Momentum.Eta() != 1.2


def test_preserves_pt(run_module):
    config = make_config("AngularSmearing", EtaResolutionFormula=0.05, PhiResolutionFormula=0.05)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    smeared = output.At(0)
    assert smeared.Momentum.Pt() == pytest.approx(50.0, rel=1e-6)
