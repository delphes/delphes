from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config("LLPFilter", {"PdgCode": [1000022]}, **extra)


def test_keeps_matching_llp(run_module):
    output = run_module(make_module_config(), [{"pid": 1000022, "charge": 1, "D1": -1, "D2": -1}])
    assert output.GetEntries() == 1


def test_rejects_non_matching(run_module):
    output = run_module(make_module_config(), [{"pid": 13, "charge": 1, "D1": -1, "D2": -1}])
    assert output.GetEntries() == 0


def test_pt_min_filter(run_module):
    output = run_module(
        make_module_config(PTMin=100.0), [{"pid": 1000022, "charge": 1, "D1": -1, "D2": -1, "pt": 50.0}]
    )
    assert output.GetEntries() == 0


def test_daughter_number_filter(run_module):
    output = run_module(make_module_config(DaughterNumber=2), [{"pid": 1000022, "charge": 1, "D1": -1, "D2": -1}])
    assert output.GetEntries() == 0


def test_empty_input(run_module):
    output = run_module(make_module_config(), [])
    assert output.GetEntries() == 0


def test_pt_at_min_kept(run_module):
    output = run_module(
        make_module_config(PTMin=100.0),
        [{"pid": 1000022, "charge": 1, "D1": -1, "D2": -1, "pt": 100.0}],
    )
    assert output.GetEntries() == 1


def test_daughter_number_at_cut_kept(run_module):
    output = run_module(
        make_module_config(DaughterNumber=2),
        [{"pid": 1000022, "charge": 1, "D1": 0, "D2": 1}],
    )
    assert output.GetEntries() == 1


def test_daughter_number_below_cut_dropped(run_module):
    output = run_module(
        make_module_config(DaughterNumber=2),
        [{"pid": 1000022, "charge": 1, "D1": 0, "D2": 0}],
    )
    assert output.GetEntries() == 0


def decay_region_config(**extra):
    return build_config(
        "LLPFilter",
        {
            "PdgCode": [1000022],
            "RequireDecayRegion": True,
            "DecayRegionEtaMax": 1.0,
            "DecayRegionEtaMin": 0.1,
            "DecayRegionZMax": 100.0,
            "DecayRegionZMin": 10.0,
            "DecayRegionRMax": 100.0,
            "DecayRegionRMin": 10.0,
        },
        **extra,
    )


def run_decay_region_test(run_generic, config):
    def setup(module, factory):
        arr = module.ExportArray("inputParticles")
        c = make_candidate(factory, 100.0, 0.5, pid=1000022, charge=1)
        c.DecayPosition.SetXYZT(30.0, 40.0, 50.0, 0.0)
        arr.Add(c)

    return run_generic(config, setup=setup)


def test_decay_region_interior_kept(run_generic):
    assert run_decay_region_test(run_generic, decay_region_config()).GetEntries() == 1


def test_decay_region_zmax_strict(run_generic):
    assert run_decay_region_test(run_generic, decay_region_config(DecayRegionZMax=50.0)).GetEntries() == 0


def test_decay_region_zmin_strict(run_generic):
    assert run_decay_region_test(run_generic, decay_region_config(DecayRegionZMin=50.0)).GetEntries() == 0


def test_decay_region_rmax_strict(run_generic):
    assert run_decay_region_test(run_generic, decay_region_config(DecayRegionRMax=50.0)).GetEntries() == 0


def test_decay_region_rmin_strict(run_generic):
    assert run_decay_region_test(run_generic, decay_region_config(DecayRegionRMin=50.0)).GetEntries() == 0


def run_llp_test(run_generic, config, pt, eta, decay=None):
    def setup(module, factory):
        arr = module.ExportArray("inputParticles")
        c = make_candidate(factory, pt, eta, pid=1000022, charge=1)
        if decay is not None:
            c.DecayPosition.SetXYZT(*decay)
        arr.Add(c)

    return run_generic(config, setup=setup)


def test_extreme_eta_outside_decay_region_dropped(run_generic):
    config = decay_region_config(DecayRegionEtaMax=2.0)
    assert run_llp_test(run_generic, config, 100.0, 5.0, decay=(30.0, 40.0, 50.0, 0.0)).GetEntries() == 0
    assert run_llp_test(run_generic, config, 100.0, -5.0, decay=(30.0, 40.0, 50.0, 0.0)).GetEntries() == 0


def test_extreme_eta_without_decay_region_kept(run_generic):
    assert run_llp_test(run_generic, make_module_config(), 100.0, 5.0).GetEntries() == 1
    assert run_llp_test(run_generic, make_module_config(), 1e4, -5.0).GetEntries() == 1


def test_near_zero_pt_dropped_by_ptmin(run_generic):
    assert run_llp_test(run_generic, make_module_config(PTMin=10.0), 1e-6, 0.5).GetEntries() == 0
