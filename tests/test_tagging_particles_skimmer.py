from conftest import add_generator_particles, add_input_partons, build_config


def test_skims_partons_with_tau(run_generic):
    config = build_config(
        "TaggingParticlesSkimmer",
        {
            "PartonInputArray": "Delphes/inputPartons",
            "ParticleInputArray": "Delphes/inputParticles",
            "OutputArray": "taggingParticles",
            "PTMin": 0.0,
            "EtaMax": 10.0,
        },
    )

    def setup(module, factory):
        add_input_partons(module, factory, [{"pt": 50.0, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}])
        add_generator_particles(
            module,
            factory,
            [{"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1}],
        )

    output = run_generic(config, setup=setup, outputs=("TestModule/taggingParticles",))

    assert output.GetEntries() == 1
    assert output.At(0).PID == 15


def test_empty_input(run_generic):
    config = build_config(
        "TaggingParticlesSkimmer",
        {
            "PartonInputArray": "Delphes/inputPartons",
            "ParticleInputArray": "Delphes/inputParticles",
            "OutputArray": "taggingParticles",
            "PTMin": 0.0,
            "EtaMax": 10.0,
        },
    )

    def setup(module, factory):
        module.ExportArray("inputPartons")
        module.ExportArray("inputParticles")

    output = run_generic(config, setup=setup, outputs=("TestModule/taggingParticles",))
    assert output.GetEntries() == 0


def make_config(**extra):
    return build_config(
        "TaggingParticlesSkimmer",
        {
            "PartonInputArray": "Delphes/inputPartons",
            "ParticleInputArray": "Delphes/inputParticles",
            "OutputArray": "taggingParticles",
            "PTMin": 10.0,
            "EtaMax": 2.5,
        },
        **extra,
    )


def tau_setup(module, factory, extra_partons=()):
    partons = [{"pt": 50.0, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}]
    partons.extend(extra_partons)
    add_input_partons(module, factory, partons)
    add_generator_particles(
        module,
        factory,
        [{"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1}],
    )


def test_no_tau_early_return_drops_other_partons(run_generic):
    def setup(module, factory):
        add_input_partons(module, factory, [{"pt": 50.0, "eta": 0.5, "pid": 5, "status": 3}])
        module.ExportArray("inputParticles")

    output = run_generic(make_config(), setup=setup, outputs=("TestModule/taggingParticles",))
    assert output.GetEntries() == 0


def test_with_tau_other_partons_kept(run_generic):
    def setup(module, factory):
        tau_setup(module, factory, extra_partons=[{"pt": 50.0, "eta": 0.5, "pid": 5, "status": 3}])

    output = run_generic(make_config(), setup=setup, outputs=("TestModule/taggingParticles",))

    assert output.GetEntries() == 2


def test_pt_at_min_kept(run_generic):
    def setup(module, factory):
        tau_setup(module, factory, extra_partons=[{"pt": 10.0, "eta": 0.5, "pid": 2, "status": 3}])

    output = run_generic(make_config(), setup=setup, outputs=("TestModule/taggingParticles",))
    assert output.GetEntries() == 2


def test_eta_max_boundary(run_generic):
    def setup_kept(module, factory):
        tau_setup(module, factory, extra_partons=[{"pt": 50.0, "eta": 2.49, "pid": 2, "status": 3}])

    output = run_generic(make_config(), setup=setup_kept, outputs=("TestModule/taggingParticles",))

    assert output.GetEntries() == 2

    def setup_dropped(module, factory):
        tau_setup(module, factory, extra_partons=[{"pt": 50.0, "eta": 2.51, "pid": 2, "status": 3}])

    output = run_generic(make_config(), setup=setup_dropped, outputs=("TestModule/taggingParticles",))

    assert output.GetEntries() == 1


def test_tiny_pt_parton_dropped(run_generic):
    def setup(module, factory):
        tau_setup(module, factory, extra_partons=[{"pt": 1e-6, "eta": 0.5, "pid": 2, "status": 3}])

    output = run_generic(make_config(), setup=setup, outputs=("TestModule/taggingParticles",))
    assert output.GetEntries() == 1


def test_extreme_eta_parton_etamax_gate(run_generic):
    def setup_pos(module, factory):
        tau_setup(module, factory, extra_partons=[{"pt": 50.0, "eta": 5.0, "pid": 2, "status": 3}])

    def setup_neg(module, factory):
        tau_setup(module, factory, extra_partons=[{"pt": 50.0, "eta": -5.0, "pid": 2, "status": 3}])

    output = run_generic(make_config(), setup=setup_pos, outputs=("TestModule/taggingParticles",))
    assert output.GetEntries() == 1
    output = run_generic(make_config(), setup=setup_neg, outputs=("TestModule/taggingParticles",))
    assert output.GetEntries() == 1
    output = run_generic(make_config(EtaMax=10.0), setup=setup_pos, outputs=("TestModule/taggingParticles",))
    assert output.GetEntries() == 2
    output = run_generic(make_config(EtaMax=10.0), setup=setup_neg, outputs=("TestModule/taggingParticles",))
    assert output.GetEntries() == 2
