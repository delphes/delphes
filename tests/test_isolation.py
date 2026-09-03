from .conftest import build_config, make_candidate


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


def test_pt_ratio_max_boundary(run_generic):
    passing = run_isolation_test(
        run_generic, make_module_config(PTRatioMax=0.21), candidates=[(50.0, 0.5, 0.0)], isolations=[(10.0, 0.5, 0.0)]
    )
    assert passing.GetEntries() == 1
    rejected = run_isolation_test(
        run_generic, make_module_config(PTRatioMax=0.19), candidates=[(50.0, 0.5, 0.0)], isolations=[(10.0, 0.5, 0.0)]
    )
    assert rejected.GetEntries() == 0


def test_empty_input(run_generic):
    results = run_isolation_test(run_generic, make_module_config(), candidates=[], isolations=[])
    assert results.GetEntries() == 0


def test_deltar_at_max_included(run_generic):
    at_cone = run_isolation_test(
        run_generic, make_module_config(), candidates=[(50.0, 0.0, 0.0)], isolations=[(10.0, 0.0, 0.5)]
    )
    assert at_cone.GetEntries() == 0
    below = run_isolation_test(
        run_generic, make_module_config(DeltaRMax=0.49), candidates=[(50.0, 0.0, 0.0)], isolations=[(10.0, 0.0, 0.5)]
    )
    assert below.GetEntries() == 1


def test_ratio_at_ptratio_max_passes(run_generic):
    results = run_isolation_test(
        run_generic, make_module_config(PTRatioMax=0.2), candidates=[(50.0, 0.5, 0.0)], isolations=[(10.0, 0.5, 0.0)]
    )
    assert results.GetEntries() == 1


def test_isolating_pt_at_min_included(run_generic):
    at_min = run_isolation_test(
        run_generic, make_module_config(PTMin=10.0), candidates=[(50.0, 0.5, 0.0)], isolations=[(10.0, 0.5, 0.0)]
    )
    assert at_min.GetEntries() == 0
    below = run_isolation_test(
        run_generic, make_module_config(PTMin=10.0), candidates=[(50.0, 0.5, 0.0)], isolations=[(9.9, 0.5, 0.0)]
    )
    assert below.GetEntries() == 1


def test_mini_cone_excludes_close_particles(run_generic):
    close = run_isolation_test(
        run_generic,
        make_module_config(UseMiniCone=True, DeltaRMin=0.01),
        candidates=[(50.0, 0.5, 0.0)],
        isolations=[(10.0, 0.5, 0.0)],
    )
    assert close.GetEntries() == 1
    inside = run_isolation_test(
        run_generic,
        make_module_config(UseMiniCone=True, DeltaRMin=0.01),
        candidates=[(50.0, 0.5, 0.0)],
        isolations=[(10.0, 0.5, 0.2)],
    )
    assert inside.GetEntries() == 0
