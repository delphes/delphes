from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "TruthVertexFinder",
        {
            "InputArray": "Delphes/inputParticles",
            "VertexOutputArray": "vertices",
            "Resolution": 0.001,
        },
        **extra,
    )


def run_vertex_test(run_generic, config, particles):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for x, y, z, pt, pid in particles:
            c = make_candidate(factory, pt, 0.5, pid=pid)
            c.Position.SetXYZT(x, y, z, 0.0)
            input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/vertices",))


def test_creates_vertex(run_generic):
    output = run_vertex_test(run_generic, make_module_config(), [(0.0, 0.0, 0.0, 50.0, 211)])
    assert output.GetEntries() >= 1


def test_same_vertex(run_generic):
    output = run_vertex_test(
        run_generic, make_module_config(), [(0.0, 0.0, 0.0, 50.0, 211), (0.0005, 0.0005, 0.0005, 30.0, 211)]
    )
    assert output.GetEntries() == 1
