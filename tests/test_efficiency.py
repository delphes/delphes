import pytest
from .conftest import assert_deterministic, make_config, run_repeated


def test_pass_all(run_module):
    config = make_config("Efficiency", EfficiencyFormula=1.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 1


def test_reject_all(run_module):
    config = make_config("Efficiency", EfficiencyFormula=0.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 0


def test_conditional_pt_cut(run_module):
    config = make_config("Efficiency", EfficiencyFormula="(pt > 10) * 1.0")
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}, {"pt": 5.0, "eta": 0.5}])
    assert output.GetEntries() == 1


def test_conditional_eta_cut(run_module):
    config = make_config("Efficiency", EfficiencyFormula="(abs(eta) < 1.5) * 1.0", UseMomentumVector=True)
    output = run_module(config, [{"pt": 50.0, "eta": 1.0}, {"pt": 50.0, "eta": 3.0}])
    assert output.GetEntries() == 1


def test_statistical_efficiency_ratio(load_delphes):
    config = make_config("Efficiency", EfficiencyFormula=0.3)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = factory.NewCandidate()
        c.Momentum.SetPtEtaPhiE(50.0, 0.5, 0.0, 50.0)
        c.Position.SetPtEtaPhiE(0.0, 0.5, 0.0, 0.0)
        input_array.Add(c)

    n = 400
    results = run_repeated(load_delphes, config, n, setup)
    kept = sum(r.GetEntries() for r in results)
    assert kept / n == pytest.approx(0.3, abs=0.08)


def test_deterministic_with_fixed_seed(run_module):
    config = make_config("Efficiency", EfficiencyFormula=0.3)
    candidates = [{"pt": 50.0, "eta": 0.5, "pid": 211, "charge": 1}]
    assert_deterministic(lambda: run_module(config, candidates), abs_tol=1e-12)


def test_empty_input(run_module):
    output = run_module(make_config("Efficiency", EfficiencyFormula=1.0), [])
    assert output.GetEntries() == 0
