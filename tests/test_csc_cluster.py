import pytest
from .conftest import assert_deterministic, candidate_snapshots, make_candidate, make_config, run_repeated


def run_cluster_test(run_generic, config, eta=3.5, ehad=50.0, eem=0.0):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, 50.0, eta, pid=130, charge=0)
        c.DecayPosition.SetXYZT(500.0, 0.0, 1000.0, 0.0)
        c.Ehad = ehad
        c.Eem = eem
        input_array.Add(c)

    return run_generic(config, setup=setup)


CSC_CLUSTERS = [
    ("CscClusterEfficiency", {}),
    ("CscClusterId", {"EtaCutFormula": 10.0, "EtaCutMax": 10.0}),
]
CSC_CLUSTER_IDS = [name for name, _ in CSC_CLUSTERS]


@pytest.mark.parametrize("class_name,extra", CSC_CLUSTERS, ids=CSC_CLUSTER_IDS)
class TestCscCluster:
    def test_pass_all(self, run_generic, class_name, extra):
        config = make_config(class_name, EfficiencyFormula=1.0, **extra)
        output = run_cluster_test(run_generic, config)
        assert output.GetEntries() == 1

    def test_reject_all(self, run_generic, class_name, extra):
        config = make_config(class_name, EfficiencyFormula=0.0, **{key: 0.0 for key in extra})
        output = run_cluster_test(run_generic, config)
        assert output.GetEntries() == 0


@pytest.mark.parametrize("class_name,extra", CSC_CLUSTERS, ids=CSC_CLUSTER_IDS)
def test_intermediate_efficiency(load_delphes, class_name, extra):
    if class_name == "CscClusterId":
        extra = dict(extra, EtaCutFormula=1.0)
    config = make_config(class_name, EfficiencyFormula=0.5, **extra)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, 50.0, 3.5, pid=130, charge=0)
        c.DecayPosition.SetXYZT(500.0, 0.0, 1000.0, 0.0)
        c.Ehad = 50.0
        c.Eem = 0.0
        input_array.Add(c)

    results = run_repeated(load_delphes, config, 100, setup)
    accepted = [1 for r in results if r.GetEntries() == 1]
    fraction = len(accepted) / len(results)
    assert fraction == pytest.approx(0.5, abs=0.2)


@pytest.mark.parametrize("class_name,extra", CSC_CLUSTERS, ids=CSC_CLUSTER_IDS)
def test_empty_input(run_generic, class_name, extra):
    def setup(module, factory):
        module.ExportArray("inputParticles")

    output = run_generic(make_config(class_name, EfficiencyFormula=1.0, **extra), setup=setup)
    assert output.GetEntries() == 0


@pytest.mark.parametrize("class_name,extra", CSC_CLUSTERS, ids=CSC_CLUSTER_IDS)
def test_deterministic_with_fixed_seed(run_generic, class_name, extra):
    if class_name == "CscClusterId":
        extra = dict(extra, EtaCutFormula=1.0)
    config = make_config(class_name, EfficiencyFormula=0.5, **extra)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for eta in (3.4, 3.5, 3.6):
            c = make_candidate(factory, 50.0, eta, pid=130, charge=0)
            c.DecayPosition.SetXYZT(500.0, 0.0, 1000.0, 0.0)
            c.Ehad = 50.0
            c.Eem = 0.0
            input_array.Add(c)

    assert_deterministic(
        lambda: run_generic(config, setup=setup),
        extract=lambda output: candidate_snapshots(output, ("PID", "Momentum")),
    )
