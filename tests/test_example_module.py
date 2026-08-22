from conftest import make_config


def test_pass_all(run_module):
    config = make_config("ExampleModule", EfficiencyFormula=1.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 1


def test_reject_all(run_module):
    config = make_config("ExampleModule", EfficiencyFormula=0.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 0


def test_conditional(run_module):
    config = make_config("ExampleModule", EfficiencyFormula="{(pt > 10) * 1.0}")
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 1
