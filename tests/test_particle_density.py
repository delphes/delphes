from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "ParticleDensity",
        {
            "InputArray": "Delphes/inputParticles",
            "OutputArray": "outputParticles",
            "UseMomentumVector": False,
            "EtaBins": [-5.0, -3.0, -1.0, 0.0, 1.0, 3.0, 5.0],
            "PhiBins": [-3.14159, -1.5708, 0.0, 1.5708, 3.14159],
        },
        **extra,
    )


def run_density_test(run_generic, config, particles):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for pt, eta, phi in particles:
            c = make_candidate(factory, pt, eta, phi, pid=22)
            c.Position.SetPtEtaPhiE(1000.0, eta, phi, pt)
            input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/outputParticles",))


def test_sets_density(run_generic):
    output = run_density_test(run_generic, make_module_config(), [(10.0, 0.5, 0.0), (10.0, 0.5, 0.1)])
    assert output.GetEntries() == 2
    assert output.At(0).ParticleDensity > 0


def test_single_particle_density(run_generic):
    output = run_density_test(run_generic, make_module_config(), [(10.0, 0.5, 0.0)])
    assert output.GetEntries() == 1


def test_different_eta_bins(run_generic):
    config = make_module_config(EtaBins=[-2.5, 0.0, 2.5], PhiBins=[-3.14159, 0.0, 3.14159])
    output = run_density_test(run_generic, config, [(10.0, 0.5, 0.0)])
    assert output.GetEntries() == 1
