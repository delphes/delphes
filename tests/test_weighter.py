from .conftest import make_config, make_particle


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


def test_empty_input(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[13, -13], 3.0])
    output = run_weighter_test(run_generic, config, [])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 1.0


def test_zero_weight(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[13], 0.0])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 13, 3)])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 0.0
    assert output.At(0).Momentum.E() == 0.0


def test_negative_weight(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[13], -2.5])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 13, 3)])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.E() == -2.5
    assert output.At(0).Momentum.Px() == 2.5
    assert output.At(0).Momentum.Pt() == 2.5


def test_code_order_independent(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[211, 13], 3.0])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 13, 3), (40.0, 0.6, 211, 3)])
    assert output.At(0).Momentum.Pt() == 3.0


def test_three_code_entry(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[13, 211, 2212], 4.0])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 13, 3), (40.0, 0.6, 211, 3), (30.0, 0.7, 2212, 3)])
    assert output.At(0).Momentum.Pt() == 4.0


def test_partial_code_set_falls_back_to_default(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[13, 211, 2212], 4.0])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 13, 3), (40.0, 0.6, 211, 3)])
    assert output.At(0).Momentum.Pt() == 1.0


def test_five_distinct_codes_fall_back_to_default(run_generic):
    config = make_config(
        "Weighter",
        OutputArray="weight",
        Weight=[[13], 2.0, [211], 2.0, [2212], 2.0, [311], 2.0, [321], 2.0],
    )
    output = run_weighter_test(
        run_generic,
        config,
        [(50.0, 0.5, 13, 3), (40.0, 0.6, 211, 3), (30.0, 0.7, 2212, 3), (20.0, 0.8, 311, 3), (10.0, 0.9, 321, 3)],
    )
    assert output.At(0).Momentum.Pt() == 1.0


def test_non_status3_particles_ignored(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[13], 2.0])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 13, 1)])
    assert output.At(0).Momentum.Pt() == 1.0


def test_status3_particle_outside_weight_set_ignored(run_generic):
    config = make_config("Weighter", OutputArray="weight", Weight=[[13], 2.0])
    output = run_weighter_test(run_generic, config, [(50.0, 0.5, 2212, 3)])
    assert output.At(0).Momentum.Pt() == 1.0
