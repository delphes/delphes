import pytest
from conftest import make_candidate, make_config


def make_module_config(**extra):
    return make_config("EnergySmearing", **{"ResolutionFormula": 0.1, **extra})


def run_module_test(run_generic, config, pt=50.0, eta=0.5, energy=50.0):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, pt, eta)
        c.Position.SetPtEtaPhiE(pt, eta, 0.0, energy)
        input_array.Add(c)

    return run_generic(config, setup=setup)


def test_zero_resolution(run_module):
    config = make_config("EnergySmearing", ResolutionFormula=0.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 1
    smeared = output.At(0)
    assert smeared.Momentum.Pt() == pytest.approx(50.0, rel=1e-6)


def test_sets_track_resolution(run_generic):
    output = run_module_test(run_generic, make_module_config())
    smeared = output.At(0)
    assert smeared.TrackResolution > 0
    assert smeared.TrackResolution < 1.0


def test_smears_energy(run_generic):
    output = run_module_test(run_generic, make_module_config())
    smeared = output.At(0)
    assert smeared.Momentum.Pt() != 50.0
