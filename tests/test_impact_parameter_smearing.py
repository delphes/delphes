from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "ImpactParameterSmearing",
        {
            "InputArray": "Delphes/inputTracks",
            "OutputArray": "tracks",
            "ResolutionFormula": 0.001,
        },
        **extra,
    )


def run_module_test(run_generic, config):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        c = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
        c.D0 = 0.0
        c.ErrorD0 = 0.001
        c.DZ = 0.0
        c.TrackResolution = 0.01

        mother = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
        c.AddCandidate(mother)

        input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/tracks",))


def test_smears_impact_parameters(run_generic):
    output = run_module_test(run_generic, make_module_config())
    assert output.GetEntries() == 1
    smeared = output.At(0)
    assert smeared.D0 != 0.0


def test_zero_resolution_preserves(run_generic):
    output = run_module_test(run_generic, make_module_config(ResolutionFormula=0.0))
    smeared = output.At(0)
    assert smeared.D0 == 0.0


def test_sets_error_d0(run_generic):
    output = run_module_test(run_generic, make_module_config(ResolutionFormula=0.002))
    smeared = output.At(0)
    assert smeared.ErrorD0 != 0.001
