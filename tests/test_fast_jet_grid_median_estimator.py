import pytest
from conftest import make_candidate, make_config


def run_rho_test(run_generic, config, particles):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for pt, eta, phi in particles:
            c = make_candidate(factory, pt, eta, phi, pid=22, charge=0)
            input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/rho",))


def make_config_grid(**extra):
    return make_config(
        "FastJetGridMedianEstimator",
        InputArray="Delphes/inputParticles",
        RhoOutputArray="rho",
        GridRange=[-5.0, 5.0, 1.0, 1.0],
        **extra,
    )


def test_computes_rho(run_generic):
    output = run_rho_test(run_generic, make_config_grid(), [(10.0, 0.5, 0.0), (10.0, 1.0, 0.0)])
    assert output.GetEntries() == 1


def test_empty_input(run_generic):
    output = run_rho_test(run_generic, make_config_grid(), [])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.E() == pytest.approx(0.0, abs=1e-9)


def test_uniform_grid_rho_value(run_generic):
    particles = []
    for eta in (-0.5, 0.5):
        for phi in (-2.5, -1.5, -0.5, 0.5, 1.5, 2.5):
            particles.append((10.0, eta, phi))
    output = run_rho_test(
        run_generic,
        make_config(
            "FastJetGridMedianEstimator",
            InputArray="Delphes/inputParticles",
            RhoOutputArray="rho",
            GridRange=[-1.0, 1.0, 1.0, 1.0],
        ),
        particles,
    )
    assert output.GetEntries() == 1
    rho = output.At(0)
    assert 9.0 < rho.Momentum.E() <= 10.0
    assert rho.Edges[0] == pytest.approx(-1.0, abs=1e-6)
    assert rho.Edges[1] == pytest.approx(1.0, abs=1e-6)
