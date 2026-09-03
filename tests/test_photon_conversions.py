import pytest
from .conftest import assert_deterministic, build_config, candidate_snapshots, make_candidate, run_repeated


def make_module_config(**extra):
    return build_config(
        "PhotonConversions",
        {
            "Radius": 1.0,
            "HalfLength": 3.0,
            "EtaMin": 2.0,
            "EtaMax": 5.0,
            "Step": 0.1,
            "ConversionMap": 0.0,
        },
        **extra,
    )


def run_conversion_test(run_generic, config, pid=22, pt=50.0, eta=3.0):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, pt, eta, pid=pid)
        input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/outputParticles",))


def test_non_photon_passes_through(run_generic):
    output = run_conversion_test(run_generic, make_module_config(), pid=211, pt=50.0, eta=3.0)
    assert output.GetEntries() == 1
    assert output.At(0).PID == 211


def test_photon_with_zero_conversion(run_generic):
    output = run_conversion_test(run_generic, make_module_config(), pid=22, pt=50.0, eta=3.0)
    assert output.GetEntries() == 1
    assert output.At(0).PID == 22


def test_photon_below_eta_min_dropped(run_generic):
    output = run_conversion_test(run_generic, make_module_config(), pid=22, pt=50.0, eta=1.5)
    assert output.GetEntries() == 0


def test_photon_above_eta_max_dropped(run_generic):
    output = run_conversion_test(run_generic, make_module_config(), pid=22, pt=50.0, eta=5.5)
    assert output.GetEntries() == 0


def test_photon_outside_cylinder_dropped(run_generic):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, 50.0, 3.0, pid=22)

        c.Position.SetXYZT(2000.0, 0.0, 0.0, 0.0)
        input_array.Add(c)

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputParticles",))
    assert output.GetEntries() == 0


def test_high_conversion_rate_produces_pair(run_generic):
    output = run_conversion_test(run_generic, make_module_config(ConversionMap=100.0), pid=22, pt=50.0, eta=3.0)
    assert output.GetEntries() == 2
    pids = sorted(output.At(i).PID for i in range(2))
    assert pids == [-11, 11]
    assert all(output.At(i).IsFromConversion == 1 for i in range(2))


def test_intermediate_conversion_probability(load_delphes):
    config = make_module_config(ConversionMap=0.3)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        input_array.Add(make_candidate(factory, 50.0, 3.0, pid=22))

    n = 400
    results = run_repeated(load_delphes, config, n, setup)
    converted = sum(1 for r in results if r.GetEntries() == 2)
    fraction = converted / n

    assert 0.2 < fraction < 0.8
    assert fraction == pytest.approx(0.5, abs=0.08)


def test_converted_pair_conserves_photon_four_momentum(load_delphes):
    config = make_module_config(ConversionMap=100.0)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, 50.0, 3.0, pid=22)
        input_array.Add(c)
        return (c.Momentum.Px(), c.Momentum.Py(), c.Momentum.Pz(), c.Momentum.E())

    module, factory, _ = load_delphes(config)
    parent = setup(module, factory)
    module.Init()
    module.Process()
    out = module.ImportArray("TestModule/outputParticles")
    assert out.GetEntries() == 2
    s_px = out.At(0).Momentum.Px() + out.At(1).Momentum.Px()
    s_py = out.At(0).Momentum.Py() + out.At(1).Momentum.Py()
    s_pz = out.At(0).Momentum.Pz() + out.At(1).Momentum.Pz()
    s_e = out.At(0).Momentum.E() + out.At(1).Momentum.E()
    assert s_px == pytest.approx(parent[0], rel=1e-9)
    assert s_py == pytest.approx(parent[1], rel=1e-9)
    assert s_pz == pytest.approx(parent[2], rel=1e-9)
    assert s_e == pytest.approx(parent[3], rel=1e-9)


def test_deterministic_with_fixed_seed(run_generic):
    config = make_module_config(ConversionMap=100.0)
    assert_deterministic(
        lambda: run_conversion_test(run_generic, config, pid=22, pt=50.0, eta=3.0),
        extract=lambda out: candidate_snapshots(out, ("PID", "Momentum", "Position")),
        abs_tol=1e-12,
    )


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputParticles")

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputParticles",))
    assert output.GetEntries() == 0
