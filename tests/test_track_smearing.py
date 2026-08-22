import pytest
from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "TrackSmearing",
        {
            "InputArray": "Delphes/inputTracks",
            "OutputArray": "outputParticles",
            "Bz": 3.8,
            "D0ResolutionFormula": "{0.001}",
            "DZResolutionFormula": "{0.002}",
            "PResolutionFormula": "{0.01}",
            "CtgThetaResolutionFormula": "{0.003}",
            "PhiResolutionFormula": "{0.004}",
        },
        **extra,
    )


def run_module_test(run_generic, config, pt=50.0, eta=0.5):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        c = make_candidate(factory, pt, eta, pid=211, charge=1)
        c.D0 = 0.0
        c.DZ = 0.0
        c.P = pt
        c.CtgTheta = 1.0
        c.Phi = 0.0
        c.IsPU = 0
        c.InitialPosition.SetXYZT(0.0, 0.0, 0.0, 0.0)
        input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/outputParticles",))


def test_smears_track(run_generic):
    output = run_module_test(run_generic, make_module_config(), pt=50.0, eta=0.5)
    assert output.GetEntries() == 1


def test_sets_errors(run_generic):
    output = run_module_test(run_generic, make_module_config(), pt=50.0, eta=0.5)
    smeared = output.At(0)
    assert smeared.ErrorD0 == pytest.approx(0.001, rel=1e-3)
    assert smeared.ErrorDZ == pytest.approx(0.002, rel=1e-3)
    assert smeared.ErrorCtgTheta == pytest.approx(0.003, rel=1e-3)
    assert smeared.ErrorPhi == pytest.approx(0.004, rel=1e-3)


def test_sets_track_resolution(run_generic):
    output = run_module_test(run_generic, make_module_config(), pt=50.0, eta=0.5)
    smeared = output.At(0)
    assert smeared.TrackResolution > 0
