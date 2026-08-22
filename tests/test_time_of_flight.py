import pytest
from conftest import build_config, make_candidate, make_vertex


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

        mother = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
        mother.Position.SetXYZT(1000.0, 0.0, 0.0, vertex_time * 299.792458)
        c.AddCandidate(mother)

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
