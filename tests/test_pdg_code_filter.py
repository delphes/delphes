from conftest import build_config


def make_module_config(**extra):
    return build_config("PdgCodeFilter", {"PdgCode": ["13"]}, **extra)


def test_removes_matching_pdg(run_module):
    output = run_module(make_module_config(), [{"pid": 13, "charge": 1}])
    assert output.GetEntries() == 0


def test_keeps_non_matching_pdg(run_module):
    output = run_module(make_module_config(), [{"pid": 11, "charge": 1}])
    assert output.GetEntries() == 1


def test_pt_min_filter(run_module):
    output = run_module(make_module_config(PTMin=100.0), [{"pid": 11, "charge": 1, "pt": 50.0}])
    assert output.GetEntries() == 0


def test_pt_min_pass(run_module):
    output = run_module(make_module_config(PTMin=10.0), [{"pid": 11, "charge": 1, "pt": 50.0}])
    assert output.GetEntries() == 1


def test_invert_keeps_matching(run_module):
    output = run_module(make_module_config(Invert=True), [{"pid": 13, "charge": 1}])
    assert output.GetEntries() == 1


def test_require_status(run_module):
    output = run_module(make_module_config(RequireStatus=True, Status=3), [{"pid": 13, "charge": -1, "status": 1}])
    assert output.GetEntries() == 0


def test_require_charge(run_module):
    output = run_module(make_module_config(RequireCharge=True, Charge=-1), [{"pid": 13, "charge": 1}])
    assert output.GetEntries() == 0


def test_require_not_pileup(run_module):
    output = run_module(make_module_config(RequireNotPileup=True), [{"pid": 13, "charge": -1, "IsPU": 1}])
    assert output.GetEntries() == 0


def test_multiple_pdg_codes(run_module):
    config = make_module_config(PdgCode=["11", "13"])
    output = run_module(
        config,
        [
            {"pt": 50.0, "eta": 0.5, "pid": 11, "charge": -1},
            {"pt": 30.0, "eta": 0.5, "pid": 211, "charge": 1},
        ],
    )
    assert output.GetEntries() == 1
