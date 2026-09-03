import math

from .conftest import make_candidate, make_config


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


def test_empty_input(run_module):
    output = run_module(make_config("ExampleModule", EfficiencyFormula=1.0), [])
    assert output.GetEntries() == 0


def test_efficiency_formula_of_pt_zero_dropped(run_module):
    config = make_config("ExampleModule", EfficiencyFormula="pt")
    output = run_module(config, [{"pt": 0.0, "eta": 0.5}])
    assert output.GetEntries() == 0


def test_efficiency_formula_of_pt_above_one_kept(run_module):
    config = make_config("ExampleModule", EfficiencyFormula="pt")
    output = run_module(config, [{"pt": 1e4, "eta": 5.0}])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 1e4


def test_formula_branches_on_position_eta(load_delphes):
    config = make_config("ExampleModule", EfficiencyFormula="{(eta > 1.0) * 1.0}")

    def run(position_eta):
        module, factory, _ = load_delphes(config)
        array = module.ExportArray("inputParticles")
        c = make_candidate(factory, 50.0, 2.0)
        c.Position.SetPtEtaPhiE(1.0, position_eta, 0.0, math.cosh(position_eta))
        array.Add(c)
        module.Init()
        module.Process()
        return module.ImportArray("TestModule/outputParticles")

    assert run(2.0).GetEntries() == 1
    assert run(0.2).GetEntries() == 0
