import math

import ROOT
from .conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "TrackCovariance",
        {
            "InputArray": "Delphes/inputTracks",
            "OutputArray": "tracks",
            "Bz": 3.8,
            "NMinHits": 6,
            "DetectorGeometry": """
                # name   z_min z_max r    th     X0   nm  st_up st_dn   res_up res_dn flag
                1 PIPE  -100   100   0.01 0.002  0.5  0   0     0       0      0      0
                1 VTXLO -0.1   0.1   0.02 0.0003 0.1  2   0     1.5708  1e-5   1e-5   1
                1 VTXLO -0.2   0.2   0.03 0.0003 0.1  2   0     1.5708  1e-5   1e-5   1
                1 VTXLO -0.3   0.3   0.04 0.0003 0.1  2   0     1.5708  1e-5   1e-5   1
                1 VTXHI -0.2   0.2   0.1  0.0004 0.1  2   0     1.5708  1e-5   1e-5   1
                1 VTXHI -0.4   0.4   0.2  0.0004 0.1  2   0     1.5708  1e-5   1e-5   1
                1 DCH   -2.0   2.0   0.3  0.02   500  1   0.01  0       1e-4   0      1
                1 DCH   -2.0   2.0   0.4  0.02   500  1  -0.02  0       1e-4   0      1
                1 DCH   -2.0   2.0   0.5  0.02   500  1   0.03  0       1e-4   0      1
                1 DCH   -2.0   2.0   0.6  0.02   500  1  -0.04  0       1e-4   0      1
                1 DCH   -2.0   2.0   0.7  0.02   500  1   0.05  0       1e-4   0      1
                1 DCH   -2.0   2.0   0.8  0.02   500  1  -0.06  0       1e-4   0      1
            """,
            "ElectronScaleFactor": 1.0,
            "MuonScaleFactor": 1.0,
            "ChargedHadronScaleFactor": 1.0,
        },
        **extra,
    )


def run_track_cov_test(run_generic, config, pt=50.0, eta=0.5, pid=211, charge=1):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        c = make_candidate(factory, pt, eta, pid=pid, charge=charge)
        mother = make_candidate(factory, pt, eta, pid=pid, charge=charge)
        c.AddCandidate(mother)
        input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/tracks",))


def test_initialization_with_geometry(load_delphes):
    config = make_module_config()
    module, factory, _ = load_delphes(config)
    input_array = module.ExportArray("inputTracks")
    c = make_candidate(factory, 50.0, 0.5, pid=211, charge=1)
    input_array.Add(c)
    module.Init()


def test_process_produces_output(run_generic):
    output = run_track_cov_test(run_generic, make_module_config())
    assert output.GetEntries() == 1


def test_preserves_charge(run_generic):
    output = run_track_cov_test(run_generic, make_module_config(), charge=1)
    assert output.At(0).Charge == 1


def test_charge_determined_by_kalman(run_generic):
    output = run_track_cov_test(run_generic, make_module_config(), charge=-1, pid=-211)
    assert abs(output.At(0).Charge) == 1


def test_preserves_pid(run_generic):
    output = run_track_cov_test(run_generic, make_module_config(), pid=211)
    assert output.At(0).PID == 211


def test_errors_are_positive(run_generic):
    output = run_track_cov_test(run_generic, make_module_config())
    out = output.At(0)
    assert out.ErrorD0 > 0
    assert out.ErrorDZ > 0
    assert out.ErrorPT > 0
    assert out.ErrorPhi > 0
    assert out.ErrorCtgTheta > 0


def test_track_resolution_is_set(run_generic):
    output = run_track_cov_test(run_generic, make_module_config())
    assert output.At(0).TrackResolution > 0


def test_smeared_momentum_differs(run_generic):
    output = run_track_cov_test(run_generic, make_module_config(), pt=50.0)
    smeared_pt = output.At(0).Momentum.Pt()
    assert smeared_pt != 50.0


def test_closest_approach_set(run_generic):
    output = run_track_cov_test(run_generic, make_module_config(), pt=50.0, eta=0.5)
    out = output.At(0)
    assert out.Xd != 0.0 or out.Yd != 0.0 or out.Zd != 0.0


def test_track_covariance_matrix_set(run_generic):
    output = run_track_cov_test(run_generic, make_module_config(), pt=50.0, eta=0.5)
    out = output.At(0)
    assert out.TrackCovariance is not None


def test_muon_with_scale_factor(run_generic):
    output = run_track_cov_test(run_generic, make_module_config(), pid=13, charge=-1)
    assert output.GetEntries() == 1
    assert output.At(0).PID == 13


def test_electron_with_scale_factor(run_generic):
    output = run_track_cov_test(run_generic, make_module_config(), pid=11, charge=-1)
    assert output.GetEntries() == 1
    assert output.At(0).PID == 11


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputTracks")

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/tracks",))
    assert output.GetEntries() == 0


def cov_elements(track):
    cov = track.TrackCovariance
    n = cov.GetNrows()
    return n, [[cov(i, j) for j in range(n)] for i in range(n)]


def test_covariance_symmetric_positive_semidefinite(run_generic):
    output = run_track_cov_test(run_generic, make_module_config())
    n, rows = cov_elements(output.At(0))
    assert n == 5
    for i in range(n):
        for j in range(n):
            assert math.isfinite(rows[i][j])
            assert rows[i][j] == rows[j][i]

    m = ROOT.TMatrixDSym(n)
    for i in range(n):
        for j in range(i, n):
            m[i][j] = rows[i][j]
    eig = ROOT.TMatrixDSymEigen(m).GetEigenValues()
    values = [eig(i) for i in range(n)]
    assert min(values) >= -1e-18
    assert max(values) > 0.0

    assert all(rows[i][i] > 0.0 for i in range(n))


def test_impact_parameter_errors_shrink_with_pt(run_generic):
    low = run_track_cov_test(run_generic, make_module_config(), pt=5.0).At(0)
    high = run_track_cov_test(run_generic, make_module_config(), pt=50.0).At(0)
    assert low.ErrorD0 > high.ErrorD0 > 0.0
    assert low.ErrorDZ > high.ErrorDZ > 0.0


def test_error_pt_grows_with_pt(run_generic):
    low = run_track_cov_test(run_generic, make_module_config(), pt=5.0).At(0)
    high = run_track_cov_test(run_generic, make_module_config(), pt=50.0).At(0)
    assert high.ErrorPT > low.ErrorPT > 0.0
    ratio = high.ErrorPT / low.ErrorPT
    assert 20.0 < ratio < 200.0

    assert high.ErrorPT / high.Momentum.Pt() > low.ErrorPT / low.Momentum.Pt()


def test_low_pt_track_gives_finite_matrix(run_generic):
    out = run_track_cov_test(run_generic, make_module_config(), pt=5.0)
    assert out.GetEntries() == 1
    track = out.At(0)
    _, rows = cov_elements(track)
    assert all(math.isfinite(v) for row in rows for v in row)
    assert all(math.isfinite(x) for x in (track.ErrorD0, track.ErrorDZ, track.ErrorPT, track.ErrorPhi))
    assert track.TrackResolution > 0.0
    assert math.isfinite(track.TrackResolution)
