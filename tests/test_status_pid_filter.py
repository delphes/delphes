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


def test_empty_input(run_module):
    output = run_module(make_module_config(), [])
    assert output.GetEntries() == 0


def test_susy_pid_range_inclusive_edges(run_module):
    assert run_module(make_module_config(), [{"pid": 1000001, "charge": 1, "status": 1}]).GetEntries() == 1
    assert run_module(make_module_config(), [{"pid": 1000039, "charge": 1, "status": 1}]).GetEntries() == 1
    assert run_module(make_module_config(), [{"pid": 1000040, "charge": 1, "status": 1}]).GetEntries() == 0


def test_pythia8_status_window_strict_edges(run_module):
    assert run_module(make_module_config(), [{"pid": 1, "charge": 1, "status": 21}]).GetEntries() == 1
    assert run_module(make_module_config(), [{"pid": 1, "charge": 1, "status": 29}]).GetEntries() == 1
    assert run_module(make_module_config(), [{"pid": 1, "charge": 1, "status": 20}]).GetEntries() == 0
    assert run_module(make_module_config(), [{"pid": 1, "charge": 1, "status": 30}]).GetEntries() == 0


def test_lepton_pid_range_strict_edges(run_module):
    assert run_module(make_module_config(), [{"pid": 11, "charge": 1, "status": 1}]).GetEntries() == 1
    assert run_module(make_module_config(), [{"pid": 16, "charge": 1, "status": 1}]).GetEntries() == 1
    assert run_module(make_module_config(), [{"pid": 10, "charge": 1, "status": 1}]).GetEntries() == 0
    assert run_module(make_module_config(), [{"pid": 17, "charge": 1, "status": 1}]).GetEntries() == 0


def test_gauge_boson_range_strict_edges(run_module):
    assert run_module(make_module_config(), [{"pid": 42, "charge": 1, "status": 1}]).GetEntries() == 1
    assert run_module(make_module_config(), [{"pid": 43, "charge": 1, "status": 1}]).GetEntries() == 0


def test_photon_requires_stable_status(run_module):
    assert run_module(make_module_config(), [{"pid": 22, "charge": 1, "status": 1}]).GetEntries() == 1
    assert run_module(make_module_config(), [{"pid": 22, "charge": 1, "status": 2}]).GetEntries() == 0


def test_pt_at_min_kept(run_module):
    output = run_module(make_module_config(PTMin=100.0), [{"pid": 13, "charge": 1, "status": 3, "pt": 100.0}])
    assert output.GetEntries() == 1


def test_tiny_pt_lepton_dropped_by_ptmin(run_module):
    output = run_module(make_module_config(PTMin=10.0), [{"pid": 11, "charge": -1, "status": 1, "pt": 1e-6}])
    assert output.GetEntries() == 0


def test_tiny_pt_bquark_exempt_from_ptmin(run_module):
    output = run_module(make_module_config(PTMin=10.0), [{"pid": 5, "charge": 1, "status": 3, "pt": 1e-6}])
    assert output.GetEntries() == 1


def test_extreme_eta_high_pt_kept(run_module):
    output = run_module(make_module_config(), [{"pid": 11, "charge": -1, "status": 1, "pt": 1e4, "eta": 5.0}])
    assert output.GetEntries() == 1
    output = run_module(make_module_config(), [{"pid": 211, "charge": 1, "status": 1, "pt": 1e4, "eta": -5.0}])
    assert output.GetEntries() == 0
