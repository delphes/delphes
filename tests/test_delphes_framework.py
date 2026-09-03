import pytest
from .conftest import make_candidate, make_config


def run_framework(load_delphes, config, setup=None):
    module, factory, _ = load_delphes(config)
    result = setup(module, factory) if setup else None
    module.Init()
    module.Process()
    return module, factory, result


def add_single_particle(module, factory, pt=50.0):
    array = module.ExportArray("inputParticles")
    c = make_candidate(factory, pt, 0.5)
    c.PID = 211
    c.Charge = 1
    array.Add(c)
    return array


def test_chained_modules_output_feeds_input(load_delphes):
    config = {
        "ExecutionPath": ["EnergyScale", "Cloner"],
        "EnergyScale": {
            "Class": "EnergyScale",
            "InputArray": "Delphes/inputParticles",
            "OutputArray": "scaled",
            "ScaleFormula": 0.5,
        },
        "Cloner": {"Class": "Cloner", "InputArray": "EnergyScale/scaled", "OutputArray": "clones"},
    }

    def setup(module, factory):
        add_single_particle(module, factory, pt=50.0)

    module, _, _ = run_framework(load_delphes, config, setup=setup)
    scaled = module.ImportArray("EnergyScale/scaled")
    clones = module.ImportArray("Cloner/clones")
    assert scaled.GetEntries() == 1
    assert scaled.At(0).Momentum.Pt() == pytest.approx(25.0, rel=1e-6)

    assert clones.GetEntries() == 1
    assert clones.At(0).Momentum.Pt() == pytest.approx(25.0, rel=1e-6)


def test_three_module_chain(load_delphes):
    config = {
        "ExecutionPath": ["EnergyScale", "Cloner", "Cloner2"],
        "EnergyScale": {
            "Class": "EnergyScale",
            "InputArray": "Delphes/inputParticles",
            "OutputArray": "scaled",
            "ScaleFormula": 2.0,
        },
        "Cloner": {"Class": "Cloner", "InputArray": "EnergyScale/scaled", "OutputArray": "clones"},
        "Cloner2": {"Class": "Cloner", "InputArray": "Cloner/clones", "OutputArray": "clones2"},
    }

    def setup(module, factory):
        add_single_particle(module, factory, pt=10.0)

    module, _, _ = run_framework(load_delphes, config, setup=setup)
    clones2 = module.ImportArray("Cloner2/clones2")
    assert clones2.GetEntries() == 1
    assert clones2.At(0).Momentum.Pt() == pytest.approx(20.0, rel=1e-6)


def test_empty_execution_path(load_delphes):
    config = {"ExecutionPath": []}
    module, _, _ = run_framework(load_delphes, config)

    assert module is not None


def smeared_pts(load_delphes, seed):
    config = make_config("MomentumSmearing", ResolutionFormula="0.05")
    config["RandomSeed"] = seed

    def setup(module, factory):
        add_single_particle(module, factory, pt=50.0)

    module, _, _ = run_framework(load_delphes, config, setup=setup)
    output = module.ImportArray("TestModule/outputParticles")
    return [output.At(i).Momentum.Pt() for i in range(output.GetEntries())]


def test_random_seed_reproducibility(load_delphes):
    pts_a = smeared_pts(load_delphes, seed=42)
    pts_b = smeared_pts(load_delphes, seed=42)
    assert pts_a == pytest.approx(pts_b, rel=1e-12)


def test_different_seed_gives_different_result(load_delphes):
    pts_42 = smeared_pts(load_delphes, seed=42)
    pts_43 = smeared_pts(load_delphes, seed=43)
    assert pts_42 != pytest.approx(pts_43, rel=1e-12)
