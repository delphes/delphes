import pytest
from .conftest import assert_deterministic, build_config, candidate_snapshots, make_candidate, mean, run_repeated, stddev


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


def make_smeared_track(factory, pt, eta):
    c = make_candidate(factory, pt, eta, pid=211, charge=1)
    c.D0 = 0.0
    c.DZ = 0.0
    c.P = pt
    c.CtgTheta = 1.0
    c.Phi = 0.0
    c.IsPU = 0
    c.InitialPosition.SetXYZT(0.0, 0.0, 0.0, 0.0)
    return c


def run_module_test(run_generic, config, pt=50.0, eta=0.5):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        input_array.Add(make_smeared_track(factory, pt, eta))

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


def test_statistical_p_and_d0(load_delphes):
    config = make_module_config()

    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        input_array.Add(make_smeared_track(factory, 50.0, 0.5))

    results = run_repeated(load_delphes, config, 200, setup, output="TestModule/outputParticles")
    ps = [results[i].At(0).P for i in range(200)]
    d0s = [results[i].At(0).D0 for i in range(200)]
    assert mean(ps) == pytest.approx(50.0, abs=0.2)
    assert stddev(ps) == pytest.approx(0.5, abs=0.1)
    assert mean(d0s) == pytest.approx(0.0, abs=3.0e-4)
    assert stddev(d0s) == pytest.approx(0.001, abs=3.0e-4)


def test_statistical_dz_phi_ctgtheta(load_delphes):
    config = make_module_config()

    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        input_array.Add(make_smeared_track(factory, 50.0, 0.5))

    results = run_repeated(load_delphes, config, 200, setup, output="TestModule/outputParticles")
    dzs = [results[i].At(0).DZ for i in range(200)]
    phis = [results[i].At(0).Phi for i in range(200)]
    ctgs = [results[i].At(0).CtgTheta for i in range(200)]
    assert mean(dzs) == pytest.approx(0.0, abs=5.0e-4)
    assert stddev(dzs) == pytest.approx(0.002, abs=5.0e-4)
    assert mean(phis) == pytest.approx(0.0, abs=5.0e-3)
    assert stddev(phis) == pytest.approx(0.004, abs=5.0e-3)
    assert mean(ctgs) == pytest.approx(1.0, abs=5.0e-4)
    assert stddev(ctgs) == pytest.approx(0.003, abs=5.0e-4)


def test_deterministic_with_fixed_seed(run_generic):
    config = make_module_config()

    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        input_array.Add(make_smeared_track(factory, 50.0, 0.5))

    assert_deterministic(
        lambda: run_generic(config, setup=setup),
        extract=lambda out: candidate_snapshots(out, ("D0", "DZ", "P", "CtgTheta", "Phi", "Momentum")),
        abs_tol=1e-12,
    )


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputTracks")

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputParticles",))
    assert output.GetEntries() == 0
