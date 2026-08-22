import pytest
from conftest import make_candidate, make_config


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
