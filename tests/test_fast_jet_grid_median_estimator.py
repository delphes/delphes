from conftest import make_candidate, make_config


def run_rho_test(run_generic, config, particles):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for pt, eta, phi in particles:
            c = make_candidate(factory, pt, eta, phi, pid=22, charge=0)
            input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/rho",))


def test_computes_rho(run_generic):
    config = make_config(
        "FastJetGridMedianEstimator",
        InputArray="Delphes/inputParticles",
        RhoOutputArray="rho",
        GridRange=[-5.0, 5.0, 1.0, 1.0],
    )
    output = run_rho_test(run_generic, config, [(10.0, 0.5, 0.0), (10.0, 1.0, 0.0)])
    assert output.GetEntries() == 1
