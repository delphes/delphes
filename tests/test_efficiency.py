from conftest import make_config


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
