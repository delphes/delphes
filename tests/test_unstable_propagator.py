import math

import pytest
from conftest import assert_deterministic, build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "UnstablePropagator",
        {
            "InputArray": "Delphes/inputParticles",
            "Radius": 1.29,
            "HalfLength": 3.0,
            "Bz": 3.8,
            "Lmin": 0.001,
        },
        **extra,
    )


def setup_unstable(module, factory):
    input_array = module.ExportArray("inputParticles")

    c = make_candidate(factory, 50.0, 0.5, pid=321, charge=1)
    c.D1 = 1
    c.D2 = 2
    input_array.Add(c)

    d1 = make_candidate(factory, 20.0, 0.5, pid=211, charge=1)
    d1.Position.SetXYZT(100.0, 0.0, 0.0, 0.0)
    d1.M1 = 0
    input_array.Add(d1)

    d2 = make_candidate(factory, 20.0, 0.5, pid=-211, charge=-1)
    d2.Position.SetXYZT(200.0, 0.0, 0.0, 0.0)
    d2.M1 = 0
    input_array.Add(d2)

    return input_array, d1, d2


def test_unstable_particle_processes(load_delphes):
    module, factory, _ = load_delphes(make_module_config())
    input_array, _, _ = setup_unstable(module, factory)
    module.Init()
    module.Process()
    assert input_array.GetEntries() == 3


def test_stable_particle_passes_through(load_delphes):
    module, factory, _ = load_delphes(make_module_config())
    input_array = module.ExportArray("inputParticles")
    c = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
    input_array.Add(c)
    module.Init()
    module.Process()
    assert input_array.GetEntries() == 1


def test_daughters_repositioned(load_delphes):
    module, factory, _ = load_delphes(make_module_config())
    _, d1, d2 = setup_unstable(module, factory)
    module.Init()
    module.Process()
    assert d1.Position.X() != 100.0 or d2.Position.X() != 200.0


def test_flight_distance_pinned(load_delphes):
    module, factory, _ = load_delphes(make_module_config())
    input_array, _, _ = setup_unstable(module, factory)
    mother = input_array.At(0)
    module.Init()
    module.Process()
    assert mother.L == pytest.approx(100.0, rel=1e-9)


def test_daughters_share_propagated_position(load_delphes):
    module, factory, _ = load_delphes(make_module_config())
    _, d1, d2 = setup_unstable(module, factory)
    module.Init()
    module.Process()

    for d, (ox, oy, oz) in ((d1, (100.0, 0.0, 0.0)), (d2, (200.0, 0.0, 0.0))):
        assert (d.Position.X(), d.Position.Y(), d.Position.Z()) != (ox, oy, oz)

    assert d1.Position.X() == pytest.approx(d2.Position.X(), abs=1e-9)
    assert d1.Position.Y() == pytest.approx(d2.Position.Y(), abs=1e-9)
    assert d1.Position.Z() == pytest.approx(d2.Position.Z(), abs=1e-9)

    r = math.sqrt(d1.Position.X() ** 2 + d1.Position.Y() ** 2 + d1.Position.Z() ** 2)
    assert r == pytest.approx(100.0, rel=0.05)


def test_momenta_unchanged(load_delphes):
    def snapshot_4m(module, factory):
        input_array, _, _ = setup_unstable(module, factory)
        return [
            (c.Momentum.Px(), c.Momentum.Py(), c.Momentum.Pz(), c.Momentum.E())
            for c in (input_array.At(i) for i in range(input_array.GetEntries()))
        ]

    module_a, factory_a, _ = load_delphes(make_module_config())
    before = snapshot_4m(module_a, factory_a)
    module_a.Init()
    module_a.Process()
    input_array = module_a.ImportArray("Delphes/inputParticles")
    after = [
        (c.Momentum.Px(), c.Momentum.Py(), c.Momentum.Pz(), c.Momentum.E())
        for c in (input_array.At(i) for i in range(input_array.GetEntries()))
    ]
    assert before == pytest.approx(after, abs=1e-12)


def test_deterministic_with_fixed_seed(load_delphes):
    def run_once():
        module, factory, _ = load_delphes(make_module_config())
        input_array, d1, d2 = setup_unstable(module, factory)
        module.Init()
        module.Process()
        return (
            input_array.At(0).L,
            d1.Position.X(),
            d1.Position.Y(),
            d1.Position.Z(),
            d2.Position.X(),
            d2.Position.Y(),
            d2.Position.Z(),
        )

    assert_deterministic(run_once, extract=lambda snap: tuple(snap), abs_tol=1e-9)


def test_empty_input(load_delphes):
    module, _, _ = load_delphes(make_module_config())
    module.ExportArray("inputParticles")
    module.Init()
    module.Process()
