import math

import pytest
from .conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "ParticleDensity",
        {
            "InputArray": "Delphes/inputParticles",
            "OutputArray": "outputParticles",
            "UseMomentumVector": False,
            "EtaBins": [-5.0, -3.0, -1.0, 0.0, 1.0, 3.0, 5.0],
            "PhiBins": [-3.14159, -1.5708, 0.0, 1.5708, 3.14159],
        },
        **extra,
    )


def run_density_test(run_generic, config, particles):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for pt, eta, phi in particles:
            c = make_candidate(factory, pt, eta, phi, pid=22)
            c.Position.SetPtEtaPhiE(1000.0, eta, phi, pt)
            input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/outputParticles",))


def test_sets_density(run_generic):
    output = run_density_test(run_generic, make_module_config(), [(10.0, 0.5, 0.0), (10.0, 0.5, 0.1)])
    assert output.GetEntries() == 2
    assert output.At(0).ParticleDensity > 0


def test_density_value(run_generic):
    output = run_density_test(run_generic, make_module_config(), [(10.0, 0.5, 0.2), (10.0, 0.5, 0.3)])
    expected = 2.0 / (1.0 * math.pi / 2.0)
    assert output.At(0).ParticleDensity == pytest.approx(expected, rel=1e-3)
    assert output.At(1).ParticleDensity == pytest.approx(expected, rel=1e-3)


def test_single_particle_density(run_generic):
    output = run_density_test(run_generic, make_module_config(), [(10.0, 0.5, 0.0)])
    assert output.GetEntries() == 1


def test_empty_input(run_generic):
    output = run_density_test(run_generic, make_module_config(), [])
    assert output.GetEntries() == 0


def test_different_eta_bins(run_generic):
    config = make_module_config(EtaBins=[-2.5, 0.0, 2.5], PhiBins=[-3.14159, 0.0, 3.14159])
    output = run_density_test(run_generic, config, [(10.0, 0.5, 0.0)])
    assert output.GetEntries() == 1


def run_density_position_test(run_generic, xyz_and_pts):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for (x, y, z), pt in xyz_and_pts:
            c = make_candidate(factory, pt, 0.5, 0.0, pid=22)
            c.Position.SetXYZT(x, y, z, 0.0)
            input_array.Add(c)

    return run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputParticles",))


def test_bin_edge_particle_goes_to_next_bin(run_generic):
    x, y = 1000.0 * math.cos(1.0), -1000.0 * math.sin(1.0)
    z_half = 1000.0 * math.tanh(0.5)
    output = run_density_position_test(
        run_generic, [((x, y, -z_half), 10.0), ((x, y, 0.0), 20.0), ((x, y, z_half), 30.0)]
    )
    assert output.GetEntries() == 3
    phi_width = 1.5708
    assert output.At(0).ParticleDensity == pytest.approx(1.0 / phi_width, rel=1e-4)
    assert output.At(1).ParticleDensity == pytest.approx(2.0 / phi_width, rel=1e-4)
    assert output.At(2).ParticleDensity == pytest.approx(2.0 / phi_width, rel=1e-4)


def test_out_of_range_particle_uses_overflow_bin(run_generic):
    z = 1000.0 * math.sinh(5.5)
    single = run_density_position_test(run_generic, [((1000.0, 0.0, z), 10.0)])
    assert single.GetEntries() == 1
    d1 = single.At(0).ParticleDensity

    pair = run_density_position_test(run_generic, [((1000.0, 0.0, z), 10.0), ((1000.0, 1.0, z), 20.0)])
    assert pair.GetEntries() == 2
    assert pair.At(0).ParticleDensity == pytest.approx(2.0 * d1, rel=1e-6)
    assert pair.At(1).ParticleDensity == pytest.approx(2.0 * d1, rel=1e-6)
