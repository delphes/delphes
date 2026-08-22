from conftest import build_config


def make_module_config(**extra):
    return build_config("StatusPidFilter", {"PTMin": 0.0}, **extra)


def test_keeps_hard_scattering(run_module):
    output = run_module(make_module_config(), [{"pid": 13, "charge": 1, "status": 3}])
    assert output.GetEntries() == 1


def test_keeps_stable_lepton(run_module):
    output = run_module(make_module_config(), [{"pid": 13, "charge": 1, "status": 1}])
    assert output.GetEntries() == 1


def test_rejects_non_physical(run_module):
    output = run_module(make_module_config(), [{"pid": 1, "charge": 1, "status": 1}])
    assert output.GetEntries() == 0


def test_pt_min_filter(run_module):
    output = run_module(make_module_config(PTMin=100.0), [{"pid": 13, "charge": 1, "status": 3, "pt": 50.0}])
    assert output.GetEntries() == 0


def test_keeps_gauge_boson(run_module):
    output = run_module(make_module_config(), [{"pid": 23, "charge": 1, "status": 3}])
    assert output.GetEntries() == 1


def test_keeps_heavy_quark(run_module):
    output = run_module(make_module_config(), [{"pid": 5, "charge": 1, "status": 3}])
    assert output.GetEntries() == 1


def test_keeps_stable_photon(run_module):
    output = run_module(make_module_config(), [{"pid": 22, "charge": 1, "status": 1}])
    assert output.GetEntries() == 1


def test_rejects_stable_hadron(run_module):
    output = run_module(make_module_config(), [{"pid": 211, "charge": 1, "status": 1}])
    assert output.GetEntries() == 0


def test_keeps_susy_particle(run_module):
    output = run_module(make_module_config(), [{"pid": 1000022, "charge": 1, "status": 1}])
    assert output.GetEntries() == 1


def test_keeps_pythia8_status(run_module):
    output = run_module(make_module_config(), [{"pid": 13, "charge": 1, "status": 23}])
    assert output.GetEntries() == 1


def test_pt_min_exempts_b_quark(run_module):
    output = run_module(make_module_config(PTMin=100.0), [{"pid": 5, "charge": 1, "status": 3, "pt": 10.0}])
    assert output.GetEntries() == 1
