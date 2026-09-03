import pytest
from .conftest import assert_deterministic, build_config, candidate_snapshots, make_candidate, mean, run_repeated, stddev


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


def make_track_with_mother(factory, pt, eta):
    c = make_candidate(factory, pt, eta, pid=211, charge=1)
    c.AddCandidate(make_candidate(factory, pt, eta, pid=211, charge=1))
    return c


def run_module_test(run_generic, config):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        c = make_track_with_mother(factory, 50.0, 0.5)
        c.D0 = 0.0
        c.ErrorD0 = 0.001
        c.DZ = 0.0
        c.TrackResolution = 0.01
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


def test_statistical_d0(load_delphes):
    config = make_module_config(ResolutionFormula=0.001)

    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        input_array.Add(make_track_with_mother(factory, 50.0, 0.5))

    results = run_repeated(load_delphes, config, 200, setup, output="TestModule/tracks")
    d0s = [results[i].At(0).D0 for i in range(200)]
    assert mean(d0s) == pytest.approx(0.0, abs=3.0e-4)
    assert stddev(d0s) == pytest.approx(0.001, abs=3.0e-4)


def test_deterministic_with_fixed_seed(run_generic):
    config = make_module_config(ResolutionFormula=0.001)

    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        input_array.Add(make_track_with_mother(factory, 50.0, 0.5))

    assert_deterministic(
        lambda: run_generic(config, setup=setup, outputs=("TestModule/tracks",)),
        extract=lambda out: candidate_snapshots(out, ("D0", "DZ", "ErrorD0", "Momentum")),
        abs_tol=1e-12,
    )


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputTracks")

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/tracks",))
    assert output.GetEntries() == 0
