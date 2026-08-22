import pytest
from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "ClusterCounting",
        {
            "InputArray": "Delphes/inputTracks",
            "OutputArray": "tracks",
            "Bz": 3.8,
            "Rmin": 0.3,
            "Rmax": 2.0,
            "Zmin": -2.1,
            "Zmax": 2.1,
            "GasOption": 0,
        },
        **extra,
    )


def run_cluster_counting_test(run_generic, config, pt=10.0, eta=0.5, pid=211, charge=1):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        c = make_candidate(factory, pt, eta, pid=pid, charge=charge)
        mother = make_candidate(factory, pt, eta, pid=pid, charge=charge)
        c.AddCandidate(mother)
        input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/tracks",))


def test_initialization(load_delphes):
    config = make_module_config()
    module, factory = load_delphes(config)
    input_array = module.ExportArray("inputTracks")
    c = make_candidate(factory, 10.0, 0.5, pid=211, charge=1)
    input_array.Add(c)
    module.Init()


def test_process_produces_output(run_generic):
    output = run_cluster_counting_test(run_generic, make_module_config())
    assert output.GetEntries() == 1


def test_preserves_momentum(run_generic):
    output = run_cluster_counting_test(run_generic, make_module_config(), pt=10.0)
    assert output.At(0).Momentum.Pt() == pytest.approx(10.0, rel=1e-3)


def test_preserves_charge(run_generic):
    output = run_cluster_counting_test(run_generic, make_module_config(), charge=1)
    assert output.At(0).Charge == 1


def test_preserves_pid(run_generic):
    output = run_cluster_counting_test(run_generic, make_module_config(), pid=211)
    assert output.At(0).PID == 211


def test_multiple_tracks(run_generic):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        for pt, charge in [(5.0, 1), (10.0, 1), (50.0, -1)]:
            c = make_candidate(factory, pt, 0.5, pid=211, charge=charge)
            mother = make_candidate(factory, pt, 0.5, pid=211, charge=charge)
            c.AddCandidate(mother)
            input_array.Add(c)

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/tracks",))
    assert output.GetEntries() == 3


def test_gas_optionhelium(run_generic):
    output = run_cluster_counting_test(run_generic, make_module_config(GasOption=0))
    assert output.GetEntries() == 1


def test_gas_option_argon(run_generic):
    output = run_cluster_counting_test(run_generic, make_module_config(GasOption=2))
    assert output.GetEntries() == 1


def test_negative_charge(run_generic):
    output = run_cluster_counting_test(run_generic, make_module_config(), charge=-1, pid=-211)
    assert output.GetEntries() == 1
    assert output.At(0).Charge == -1
