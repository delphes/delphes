import pytest
from .conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "ParticlePropagator",
        {
            "OutputArray": "stableParticles",
            "NeutralOutputArray": "neutralParticles",
            "ChargedHadronOutputArray": "chargedHadrons",
            "ElectronOutputArray": "electrons",
            "MuonOutputArray": "muons",
            "Radius": 1.29,
            "HalfLength": 3.0,
            "Bz": 0.0,
        },
        **extra,
    )


def run_propagator_test(run_generic, config, particles):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for pt, eta, phi, charge, pid in particles:
            c = make_candidate(factory, pt, eta, phi, pid=pid, charge=charge)
            input_array.Add(c)

    return run_generic(
        config,
        setup=setup,
        outputs={
            "stable": "TestModule/stableParticles",
            "neutral": "TestModule/neutralParticles",
            "charged": "TestModule/chargedHadrons",
            "electrons": "TestModule/electrons",
            "muons": "TestModule/muons",
        },
    )


def test_neutral_particle_propagates(run_generic):
    results = run_propagator_test(run_generic, make_module_config(), [(50.0, 0.5, 0.0, 0, 22)])
    assert results["stable"].GetEntries() == 1
    assert results["neutral"].GetEntries() == 1
    assert results["neutral"].At(0).L > 0


def test_charged_particle_straight_line(run_generic):
    results = run_propagator_test(run_generic, make_module_config(), [(50.0, 0.5, 0.0, 1, 211)])
    assert results["stable"].GetEntries() == 1
    assert results["charged"].GetEntries() == 1
    assert results["charged"].At(0).PID == 211
    assert results["charged"].At(0).L > 0


def test_charged_particle_magnetic_field(run_generic):
    results = run_propagator_test(run_generic, make_module_config(Bz=3.8), [(50.0, 0.5, 0.0, 1, 211)])
    assert results["stable"].GetEntries() == 1
    assert results["stable"].At(0).L > 0


def test_electron_classification(run_generic):
    results = run_propagator_test(run_generic, make_module_config(), [(50.0, 0.5, 0.0, -1, 11)])
    assert results["electrons"].GetEntries() == 1
    assert results["electrons"].At(0).PID == 11


def test_muon_classification(run_generic):
    results = run_propagator_test(run_generic, make_module_config(), [(50.0, 0.5, 0.0, -1, 13)])
    assert results["muons"].GetEntries() == 1
    assert results["muons"].At(0).PID == 13


def test_particle_outside_cylinder(run_generic):
    config = make_module_config(Radius=0.5, RadiusMax=2.0, HalfLength=0.5, HalfLengthMax=2.0)

    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
        c.Position.SetXYZT(1000.0, 0.0, 0.0, 0.0)
        input_array.Add(c)

    stable = run_generic(config, setup=setup, outputs=("TestModule/stableParticles",))
    assert stable.GetEntries() == 1
    assert stable.At(0).L == 0.0


def test_trajectory_length(run_generic):
    results = run_propagator_test(run_generic, make_module_config(), [(50.0, 0.0, 0.0, 0, 22)])
    assert results["stable"].At(0).L == pytest.approx(1290.0, rel=1e-3)


def test_multiple_particles(run_generic):
    results = run_propagator_test(
        run_generic,
        make_module_config(),
        [(50.0, 0.5, 0.0, 0, 22), (50.0, -0.5, 0.0, 1, 211), (50.0, 0.0, 1.0, -1, 13)],
    )
    assert results["stable"].GetEntries() == 3
    assert results["neutral"].GetEntries() == 1
    assert results["charged"].GetEntries() == 1
    assert results["muons"].GetEntries() == 1


def test_empty_input(run_generic):
    results = run_propagator_test(run_generic, make_module_config(), [])
    for name in ("stable", "neutral", "charged", "electrons", "muons"):
        assert results[name].GetEntries() == 0
