from conftest import build_config, make_candidate


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
