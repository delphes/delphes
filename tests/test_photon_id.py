from .conftest import assert_deterministic, build_config, candidate_snapshots, make_candidate


def make_module_config(**extra):
    return build_config(
        "PhotonID",
        {
            "InputPhotonArray": "Delphes/inputPhotons",
            "InputGenArray": "Delphes/inputGen",
            "OutputArray": "photons",
            "PromptFormula": 1.0,
            "NonPromptFormula": 1.0,
            "FakeFormula": 1.0,
            "PTMin": 0.0,
            "RelIsoMax": 0.3,
        },
        **extra,
    )


def run_photonid_test(run_generic, config, recos, gens):
    def setup(module, factory):
        reco_array = module.ExportArray("inputPhotons")
        for pt, eta, iso_var in recos:
            c = make_candidate(factory, pt, eta, pid=22)
            c.Position.SetXYZT(1000.0, eta, 0.0, pt)
            c.IsolationVar = iso_var
            reco_array.Add(c)

        gen_array = module.ExportArray("inputGen")
        for pt, eta in gens:
            c = make_candidate(factory, pt, eta, pid=22)
            gen_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/photons",))


def test_prompt_photon(run_generic):
    output = run_photonid_test(run_generic, make_module_config(), recos=[(50.0, 0.5, 0.1)], gens=[(50.0, 0.5)])
    assert output.GetEntries() == 1
    assert output.At(0).Status == 1


def test_fake_photon_no_gen_match(run_generic):
    output = run_photonid_test(run_generic, make_module_config(), recos=[(50.0, 0.5, 0.1)], gens=[])
    assert output.GetEntries() == 1
    assert output.At(0).Status == 3


def test_nonprompt_photon(run_generic):
    output = run_photonid_test(run_generic, make_module_config(), recos=[(50.0, 0.5, 0.5)], gens=[(50.0, 0.5)])
    assert output.GetEntries() == 1
    assert output.At(0).Status == 2


def test_rel_iso_max_boundary(run_generic):
    at_cut = run_photonid_test(run_generic, make_module_config(), recos=[(50.0, 0.5, 0.3)], gens=[(50.0, 0.5)])
    assert at_cut.GetEntries() == 1
    assert at_cut.At(0).Status == 2
    below = run_photonid_test(run_generic, make_module_config(), recos=[(50.0, 0.5, 0.29)], gens=[(50.0, 0.5)])
    assert below.At(0).Status == 1


def test_pt_min_filter(run_generic):
    output = run_photonid_test(
        run_generic, make_module_config(PTMin=100.0), recos=[(50.0, 0.5, 0.1)], gens=[(50.0, 0.5)]
    )
    assert output.GetEntries() == 0


def test_empty_input(run_generic):
    output = run_photonid_test(run_generic, make_module_config(), recos=[], gens=[])
    assert output.GetEntries() == 0


def test_pt_at_min_passes(run_generic):
    output = run_photonid_test(
        run_generic, make_module_config(PTMin=50.0), recos=[(50.0, 0.5, 0.1)], gens=[(50.0, 0.5)]
    )
    assert output.GetEntries() == 1
    assert output.At(0).Status == 1


def test_fake_match_at_dpt_boundary(run_generic):
    output = run_photonid_test(run_generic, make_module_config(), recos=[(100.0, 0.5, 0.1)], gens=[(150.0, 0.5)])
    assert output.GetEntries() == 1
    assert output.At(0).Status == 1


def test_fake_match_above_dpt_boundary(run_generic):
    output = run_photonid_test(run_generic, make_module_config(), recos=[(100.0, 0.5, 0.1)], gens=[(150.01, 0.5)])
    assert output.GetEntries() == 1
    assert output.At(0).Status == 3


def test_fake_match_at_deltar_boundary(run_generic):
    output = run_photonid_test(run_generic, make_module_config(), recos=[(50.0, 0.5, 0.1)], gens=[(50.0, 0.6)])
    assert output.GetEntries() == 1
    assert output.At(0).Status == 3
    output = run_photonid_test(run_generic, make_module_config(), recos=[(50.0, 0.5, 0.1)], gens=[(50.0, 0.55)])
    assert output.GetEntries() == 1
    assert output.At(0).Status == 1


def test_deterministic_with_fixed_seed(run_generic):
    config = make_module_config(PromptFormula=0.5, NonPromptFormula=0.5, FakeFormula=0.5)
    recos = [(50.0, 0.5, 0.1), (50.0, 0.6, 0.1), (50.0, 0.7, 0.1)]
    gens = [(50.0, 0.5), (50.0, 0.6), (50.0, 0.7)]
    assert_deterministic(
        lambda: run_photonid_test(run_generic, config, recos, gens),
        extract=lambda out: candidate_snapshots(out, ("Status", "Momentum")),
    )
