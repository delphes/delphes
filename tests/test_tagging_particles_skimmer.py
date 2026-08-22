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
    assert output.GetEntries() >= 1
