from conftest import make_config, make_particle


def run_weighter_test(run_generic, config, particles):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for pt, eta, pid, status in particles:
            p = make_particle(factory, pt, eta, pid, status)
            input_array.Add(p)

    return run_generic(config, setup=setup, outputs=("TestModule/weight",))


def test_default_weight(run_generic):
    config = make_config("Weighter", OutputArray="weight")
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 211, 3)])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 1.0


def test_custom_weight(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[5, -5], 2.0])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 5, 3), (50.0, 0.5, -5, 3)])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 2.0


def test_no_matching_particles(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[5, -5], 2.0])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 211, 3)])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 1.0


def test_multiple_status3_particles(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[13, -13], 3.0])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 13, 3), (50.0, 0.5, -13, 3), (30.0, 1.0, 211, 3)])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 3.0
