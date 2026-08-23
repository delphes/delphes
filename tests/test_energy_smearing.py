import pytest
from conftest import (
    assert_deterministic,
    candidate_snapshots,
    make_candidate,
    make_config,
    mean,
    run_repeated,
    stddev,
)


def make_module_config(**extra):
    return make_config("EnergySmearing", **{"ResolutionFormula": 0.1, **extra})


def make_energy_particle(factory, pt, eta, energy):
    c = make_candidate(factory, pt, eta)
    c.Position.SetPtEtaPhiE(pt, eta, 0.0, energy)
    return c


def run_module_test(run_generic, config, pt=50.0, eta=0.5, energy=50.0):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        input_array.Add(make_energy_particle(factory, pt, eta, energy))

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


def test_statistical_mean_and_rms(load_delphes):
    config = make_config("EnergySmearing", ResolutionFormula=0.1)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        input_array.Add(make_energy_particle(factory, 50.0, 0.0, 50.0))

    results = run_repeated(load_delphes, config, 200, setup)
    energies = [results[i].At(0).Momentum.E() for i in range(200)]
    assert mean(energies) == pytest.approx(50.0, abs=0.05)
    assert stddev(energies) == pytest.approx(0.1, abs=0.02)


def test_deterministic_with_fixed_seed(run_generic):
    config = make_module_config()

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        input_array.Add(make_energy_particle(factory, 50.0, 0.5, 50.0))

    assert_deterministic(
        lambda: run_generic(config, setup=setup),
        extract=lambda out: candidate_snapshots(out, ("Momentum", "TrackResolution")),
        abs_tol=1e-12,
    )


def test_empty_input(run_module):
    output = run_module(make_config("EnergySmearing", ResolutionFormula=0.1), [])
    assert output.GetEntries() == 0
