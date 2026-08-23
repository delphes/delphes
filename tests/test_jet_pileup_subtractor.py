import math

import pytest
from conftest import build_config, make_jet, make_vertex


def make_module_config(**extra):
    return build_config(
        "JetPileUpSubtractor",
        {
            "JetInputArray": "Delphes/inputJets",
            "RhoInputArray": "Delphes/inputRho",
            "OutputArray": "jets",
            "JetPTMin": 0.0,
        },
        **extra,
    )


def run_subtractor_test(run_generic, config, jet_pt=50.0, rho_pt=5.0, jet_area=0.5):
    def setup(module, factory):
        jet_array = module.ExportArray("inputJets")
        jet = make_jet(factory, jet_pt, 0.5)
        jet.Area.SetPtEtaPhiE(jet_area, 0.5, 0.0, jet_area)
        jet_array.Add(jet)

        rho_array = module.ExportArray("inputRho")
        rho = make_vertex(factory)
        rho.Momentum.SetPtEtaPhiE(rho_pt, 0.0, 0.0, rho_pt)
        rho.Edges[0] = -2.5
        rho.Edges[1] = 2.5
        rho_array.Add(rho)

    return run_generic(config, setup=setup, outputs=("TestModule/jets",))


def test_subtracts_rho(run_generic):
    output = run_subtractor_test(run_generic, make_module_config(), jet_pt=50.0, rho_pt=5.0, jet_area=2.0)
    assert output.GetEntries() == 1
    subtracted_pt = output.At(0).Momentum.Pt()
    assert subtracted_pt == pytest.approx(50.0 - 5.0 * 2.0, rel=1e-3)


def test_drops_low_pt_jet(run_generic):
    output = run_subtractor_test(run_generic, make_module_config(JetPTMin=20.0), jet_pt=15.0, rho_pt=5.0, jet_area=2.0)
    assert output.GetEntries() == 0


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputJets")
        module.ExportArray("inputRho")

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/jets",))
    assert output.GetEntries() == 0


def run_full_p4_test(run_generic, jet_pt, jet_eta, rho_pt, area_pt, rho_edges, ptmin=0.0):
    def setup(module, factory):
        jet_array = module.ExportArray("inputJets")
        jet = make_jet(factory, jet_pt, jet_eta)
        jet.Area.SetPtEtaPhiE(area_pt, 0.0, 0.0, area_pt)
        jet_array.Add(jet)

        rho_array = module.ExportArray("inputRho")
        for lo, hi in rho_edges:
            rho = make_vertex(factory)
            rho.Momentum.SetPtEtaPhiE(rho_pt, 0.0, 0.0, rho_pt)
            rho.Edges[0] = lo
            rho.Edges[1] = hi
            rho_array.Add(rho)

    return run_generic(make_module_config(JetPTMin=ptmin), setup=setup, outputs=("TestModule/jets",))


def test_full_four_vector_subtraction(run_generic):
    jet_eta = 0.5
    output = run_full_p4_test(run_generic, 50.0, jet_eta, 5.0, 2.0, [(-2.5, 2.5)])
    assert output.GetEntries() == 1
    m = output.At(0).Momentum
    assert m.Px() == pytest.approx(50.0 - 5.0 * 2.0, rel=1e-9)
    assert m.Py() == pytest.approx(0.0, abs=1e-9)

    assert m.Pz() == pytest.approx(50.0 * math.sinh(jet_eta), rel=1e-9)

    assert m.E() == pytest.approx(50.0 * math.cosh(jet_eta) - 5.0 * 2.0, rel=1e-9)
    assert m.Pt() == pytest.approx(40.0, rel=1e-9)


def test_pt_equal_rho_area_is_dropped(run_generic):
    output = run_full_p4_test(run_generic, 10.0, 0.5, 5.0, 2.0, [(-2.5, 2.5)])
    assert output.GetEntries() == 0


def test_rho_selected_by_eta_range(run_generic):
    output = run_full_p4_test(run_generic, 50.0, 2.5, 5.0, 2.0, [(-2.5, 2.5)])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == pytest.approx(50.0, rel=1e-9)

    output = run_full_p4_test(run_generic, 50.0, -2.4, 5.0, 2.0, [(-2.5, 2.5)])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == pytest.approx(40.0, rel=1e-9)

    def setup_last_wins(module, factory):
        jet_array = module.ExportArray("inputJets")
        jet = make_jet(factory, 50.0, 0.5)
        jet.Area.SetPtEtaPhiE(2.0, 0.0, 0.0, 2.0)
        jet_array.Add(jet)
        rho_array = module.ExportArray("inputRho")
        for lo, hi, pt in [(-5.0, -1.0, 1.0), (-1.0, 5.0, 5.0)]:
            rho = make_vertex(factory)
            rho.Momentum.SetPtEtaPhiE(pt, 0.0, 0.0, pt)
            rho.Edges[0] = lo
            rho.Edges[1] = hi
            rho_array.Add(rho)

    output = run_generic(make_module_config(), setup=setup_last_wins, outputs=("TestModule/jets",))
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == pytest.approx(40.0, rel=1e-9)


def test_ptmin_boundary_after_subtraction(run_generic):
    assert run_full_p4_test(run_generic, 50.0, 0.5, 5.0, 8.0, [(-2.5, 2.5)], ptmin=10.0).GetEntries() == 0

    output = run_full_p4_test(run_generic, 50.5, 0.5, 5.0, 8.0, [(-2.5, 2.5)], ptmin=10.0)
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == pytest.approx(10.5, rel=1e-9)
