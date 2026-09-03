import math

import pytest
from .conftest import build_config, make_candidate, make_vertex


def make_module_config(**extra):
    return build_config(
        "TimeOfFlight",
        {
            "InputArray": "Delphes/inputTracks",
            "VertexInputArray": "Delphes/inputVertices",
            "OutputArray": "outputTracks",
            "VertexTimeMode": 1,
        },
        **extra,
    )


def run_tof_test(run_generic, config, vertex_time=0.0):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        c = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
        c.Position.SetXYZT(1000.0, 0.0, 0.0, vertex_time * 299.792458)

        gen = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
        gen.Position.SetXYZT(0.0, 0.0, 0.0, vertex_time * 299.792458)
        c.AddCandidate(gen)

        input_array.Add(c)

        vertex_array = module.ExportArray("inputVertices")
        v = make_vertex(factory, z=0.0, t=vertex_time * 299.792458)
        vertex_array.Add(v)

    return run_generic(config, setup=setup, outputs=("TestModule/outputTracks",))


def test_sets_vertex_time(run_generic):
    output = run_tof_test(run_generic, make_module_config(), vertex_time=1.0)
    assert output.GetEntries() == 1


def test_preserves_momentum(run_generic):
    output = run_tof_test(run_generic, make_module_config(), vertex_time=1.0)
    assert output.At(0).Momentum.Pt() == pytest.approx(50.0, rel=1e-6)


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputTracks")
        module.ExportArray("inputVertices")

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputTracks",))
    assert output.GetEntries() == 0


def run_mode_test(run_generic, mode, offset_mm=10.0, px=30.0, pz=40.0, e=62.0, link_vertex=True):
    config = make_module_config(VertexTimeMode=mode)

    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        c = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
        c.Momentum.SetPxPyPzE(px, 0.0, pz, e)
        c.Position.SetXYZT(1000.0, 0.0, 0.0, 1000.0)
        c.InitialPosition.SetXYZT(offset_mm, 0.0, 0.0, 0.0)

        gen = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
        gen.Momentum.SetPxPyPzE(px, 0.0, pz, e)
        gen.Position.SetXYZT(0.0, 0.0, 0.0, 0.0)
        c.AddCandidate(gen)

        input_array.Add(c)

        vertex_array = module.ExportArray("inputVertices")
        v = make_vertex(factory, z=0.0, t=0.0)
        if link_vertex:
            v.AddCandidate(gen)
        vertex_array.Add(v)

    return run_generic(config, setup=setup, outputs=("TestModule/outputTracks",))


def test_vertex_time_mode_zero_uses_truth_time(run_generic):
    output = run_tof_test(run_generic, make_module_config(VertexTimeMode=0), vertex_time=1.0)
    assert output.At(0).InitialPosition.T() == pytest.approx(299.792458, rel=1e-12)


def test_vertex_time_mode_one_zeroes_initial_time(run_generic):
    output = run_tof_test(run_generic, make_module_config(VertexTimeMode=1), vertex_time=1.0)
    assert output.At(0).InitialPosition.T() == 0.0


def test_vertex_time_mode_two_divides_by_vertex_beta(run_generic):
    output = run_mode_test(run_generic, 2, offset_mm=10.0, px=30.0, pz=40.0, e=62.0, link_vertex=True)
    assert output.GetEntries() == 1
    assert output.At(0).InitialPosition.T() == pytest.approx(10.0 * 62.0 / 50.0, rel=1e-9)


def test_vertex_time_mode_two_no_vertex_match_uses_beta_one(run_generic):
    output = run_mode_test(run_generic, 2, offset_mm=10.0, px=30.0, pz=40.0, e=62.0, link_vertex=False)
    assert output.GetEntries() == 1
    assert output.At(0).InitialPosition.T() == pytest.approx(10.0, rel=1e-12)


def test_massive_momentum_and_final_position_preserved(run_generic):
    output = run_mode_test(run_generic, 1, offset_mm=10.0, px=30.0, pz=40.0, e=62.0, link_vertex=True)
    out = output.At(0)
    assert out.Position.X() == pytest.approx(1000.0, rel=1e-12)
    assert out.Position.T() == pytest.approx(1000.0, rel=1e-12)
    assert out.Momentum.Px() == pytest.approx(30.0, rel=1e-12)
    assert out.Momentum.Pz() == pytest.approx(40.0, rel=1e-12)
    assert out.Momentum.E() == pytest.approx(62.0, rel=1e-12)
    assert out.Momentum.M() == pytest.approx(math.sqrt(62.0**2 - 50.0**2), rel=1e-12)


def test_output_keeps_reference_to_input_track(load_delphes):
    config = make_module_config()
    module, factory, _ = load_delphes(config)
    input_array = module.ExportArray("inputTracks")
    c = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
    c.Position.SetXYZT(1000.0, 0.0, 0.0, 0.0)
    gen = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
    c.AddCandidate(gen)
    input_array.Add(c)
    module.ExportArray("inputVertices").Add(make_vertex(factory))
    module.Init()
    module.Process()
    out = module.ImportArray("TestModule/outputTracks")
    assert out.GetEntries() == 1
    subs = out.At(0).GetCandidates()
    assert subs.GetEntries() == 2
    assert subs.At(0).GetUniqueID() == gen.GetUniqueID()
    assert subs.At(1).GetUniqueID() == c.GetUniqueID()
