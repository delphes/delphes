from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "Isolation",
        {
            "CandidateInputArray": "Delphes/inputCandidates",
            "IsolationInputArray": "Delphes/inputIsolation",
            "OutputArray": "outputCandidates",
            "DeltaRMax": 0.5,
            "PTRatioMax": 0.1,
            "PTMin": 0.5,
        },
        **extra,
    )


def run_isolation_test(run_generic, config, candidates, isolations):
    def setup(module, factory):
        candidate_array = module.ExportArray("inputCandidates")
        for pt, eta, phi in candidates:
            c = make_candidate(factory, pt, eta, phi, pid=11, charge=-1)
            c.Position.SetXYZT(1000.0, eta, phi, pt)
            c.IsolationVar = 0.0
            candidate_array.Add(c)

        isolation_array = module.ExportArray("inputIsolation")
        for pt, eta, phi in isolations:
            c = make_candidate(factory, pt, eta, phi, pid=211, charge=1)
            c.Position.SetXYZT(1000.0, eta, phi, pt)
            isolation_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/outputCandidates",))


def test_isolated_candidate_passes(run_generic):
    results = run_isolation_test(run_generic, make_module_config(), candidates=[(50.0, 0.5, 0.0)], isolations=[])
    assert results.GetEntries() == 1


def test_non_isolated_candidate_rejected(run_generic):
    results = run_isolation_test(
        run_generic, make_module_config(), candidates=[(50.0, 0.5, 0.0)], isolations=[(20.0, 0.5, 0.0)]
    )
    assert results.GetEntries() == 0
