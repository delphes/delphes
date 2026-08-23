import pytest
from conftest import assert_deterministic, candidate_snapshots, make_candidate, make_config, mean, run_repeated, stddev


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


def test_statistical_mean_and_rms(load_delphes):
    config = make_config("MomentumSmearing", ResolutionFormula=0.1)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        input_array.Add(make_candidate(factory, 50.0, 0.0))

    results = run_repeated(load_delphes, config, 200, setup)
    pts = [results[i].At(0).Momentum.Pt() for i in range(200)]
    assert mean(pts) == pytest.approx(50.0, abs=1.0)
    assert stddev(pts) == pytest.approx(5.0, abs=1.0)


def test_deterministic_with_fixed_seed(run_module):
    config = make_config("MomentumSmearing", ResolutionFormula=0.1)
    candidates = [{"pt": 50.0, "eta": 0.5}]
    assert_deterministic(
        lambda: run_module(config, candidates),
        extract=lambda out: candidate_snapshots(out, ("Momentum", "TrackResolution")),
        abs_tol=1e-12,
    )


def test_empty_input(run_module):
    output = run_module(make_config("MomentumSmearing", ResolutionFormula=0.1), [])
    assert output.GetEntries() == 0
