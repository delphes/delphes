import pytest
from conftest import assert_deterministic, candidate_snapshots, make_candidate, make_config, mean, run_repeated, stddev


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


def test_statistical_mean_and_rms(load_delphes):
    config = make_config("AngularSmearing", EtaResolutionFormula=0.05, PhiResolutionFormula=0.05)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        input_array.Add(make_candidate(factory, 50.0, 1.2, phi=0.3))

    results = run_repeated(load_delphes, config, 200, setup)
    etas = [results[i].At(0).Momentum.Eta() for i in range(200)]
    phis = [results[i].At(0).Momentum.Phi() for i in range(200)]
    assert mean(etas) == pytest.approx(1.2, abs=0.01)
    assert stddev(etas) == pytest.approx(0.05, abs=0.01)
    assert mean(phis) == pytest.approx(0.3, abs=0.01)
    assert stddev(phis) == pytest.approx(0.05, abs=0.01)


def test_deterministic_with_fixed_seed(run_module):
    config = make_config("AngularSmearing", EtaResolutionFormula=0.05, PhiResolutionFormula=0.05)
    candidates = [{"pt": 50.0, "eta": 1.2, "phi": 0.3}]
    assert_deterministic(
        lambda: run_module(config, candidates),
        extract=lambda out: candidate_snapshots(out, ("Momentum",)),
        abs_tol=1e-12,
    )


def test_empty_input(run_module):
    config = make_config("AngularSmearing", EtaResolutionFormula=0.05, PhiResolutionFormula=0.05)
    output = run_module(config, [])
    assert output.GetEntries() == 0
