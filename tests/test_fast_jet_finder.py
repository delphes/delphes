import math

import pytest
from .conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "FastJetFinder",
        {
            "InputArray": "Delphes/inputParticles",
            "OutputArray": "jets",
            "JetAlgorithm": 6,
            "ParameterR": 0.5,
            "JetPTMin": 5.0,
        },
        **extra,
    )


def run_jet_test(run_generic, config, particles):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        for pt, eta, phi, charge in particles:
            c = make_candidate(factory, pt, eta, phi, pid=211, charge=charge)
            input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/jets",))


def test_clusters_close_particles(run_generic):
    particles = [(20.0, 0.0, 0.0, 1), (20.0, 0.1, 0.0, 1), (20.0, 0.2, 0.0, 1)]
    jets = run_jet_test(run_generic, make_module_config(), particles)
    assert jets.GetEntries() == 1
    assert jets.At(0).Momentum.Pt() == pytest.approx(60.0, rel=1e-3)


def test_empty_input(run_generic):
    jets = run_jet_test(run_generic, make_module_config(), [])
    assert jets.GetEntries() == 0


def test_separates_distant_particles(run_generic):
    particles = [(20.0, 0.0, 0.0, 1), (30.0, 3.0, 0.0, 1)]
    jets = run_jet_test(run_generic, make_module_config(), particles)
    assert jets.GetEntries() == 2


def test_jet_pt_min_filter(run_generic):
    particles = [(20.0, 0.0, 0.0, 1), (20.0, 0.1, 0.0, 1)]
    jets = run_jet_test(run_generic, make_module_config(JetPTMin=100.0), particles)
    assert jets.GetEntries() == 0


def test_kt_algorithm(run_generic):
    particles = [(20.0, 0.0, 0.0, 1), (20.0, 0.1, 0.0, 1), (20.0, 0.2, 0.0, 1)]
    jets = run_jet_test(run_generic, make_module_config(JetAlgorithm=4), particles)
    assert jets.GetEntries() == 1


def test_cambridge_algorithm(run_generic):
    particles = [(20.0, 0.0, 0.0, 1), (20.0, 0.1, 0.0, 1), (20.0, 0.2, 0.0, 1)]
    jets = run_jet_test(run_generic, make_module_config(JetAlgorithm=5), particles)
    assert jets.GetEntries() == 1


def test_jet_constituents(run_generic):
    particles = [(20.0, 0.0, 0.0, 1), (20.0, 0.1, 0.0, 1), (20.0, 0.2, 0.0, 1), (30.0, 3.0, 0.0, 1)]
    jets = run_jet_test(run_generic, make_module_config(), particles)
    assert jets.GetEntries() == 2
    assert jets.At(0).GetCandidates().GetEntries() == 3


def test_jet_charge(run_generic):
    particles = [(20.0, 0.0, 0.0, 1), (20.0, 0.1, 0.0, -1)]
    jets = run_jet_test(run_generic, make_module_config(), particles)
    assert jets.GetEntries() == 1
    assert jets.At(0).Charge == 0


def test_jet_four_momentum_equals_constituent_sum(run_generic):
    particles = [(20.0, 0.0, 0.0, 1), (15.0, 0.3, 0.0, -1)]
    jets = run_jet_test(run_generic, make_module_config(), particles)
    assert jets.GetEntries() == 1
    jet = jets.At(0)
    exp_px = 20.0 + 15.0 * math.cos(0.0)
    exp_py = 15.0 * math.sin(0.0)
    exp_pz = 20.0 * math.sinh(0.0) + 15.0 * math.sinh(0.3)
    exp_e = 20.0 + 15.0 * math.cosh(0.3)
    assert jet.Momentum.Px() == pytest.approx(exp_px, rel=1e-9)
    assert jet.Momentum.Py() == pytest.approx(exp_py, rel=1e-9)
    assert jet.Momentum.Pz() == pytest.approx(exp_pz, rel=1e-9)
    assert jet.Momentum.E() == pytest.approx(exp_e, rel=1e-9)

    assert jet.Charge == 0
    assert jet.NCharged == 2
    assert jet.NNeutrals == 0


def test_jet_pt_min_boundary(run_generic):
    particles_at_cut = [(5.0, 0.5, 0.0, 1)]
    assert run_jet_test(run_generic, make_module_config(), particles_at_cut).GetEntries() == 1
    particles_below = [(4.9999, 0.5, 0.0, 1)]
    assert run_jet_test(run_generic, make_module_config(), particles_below).GetEntries() == 0


def test_delta_eta_phi_pin_max_constituent_spread(run_generic):
    eta2, phi2 = 0.3, 0.3
    particles = [(20.0, 0.0, 0.0, 1), (15.0, eta2, phi2, -1)]
    jets = run_jet_test(run_generic, make_module_config(), particles)
    assert jets.GetEntries() == 1
    jet = jets.At(0)

    jet = jets.At(0)
    jet_eta = jet.Momentum.Eta()
    jet_phi = jet.Momentum.Phi()
    exp_deta = max(abs(0.0 - jet_eta), abs(eta2 - jet_eta))
    exp_dphi = max(abs(0.0 - jet_phi), abs(phi2 - jet_phi))
    assert jet.DeltaEta == pytest.approx(exp_deta, abs=1e-6)
    assert jet.DeltaPhi == pytest.approx(exp_dphi, abs=1e-6)

    assert exp_deta > 0.1


def test_nsubjettiness_values_and_hierarchy(run_generic):
    config = make_module_config(ComputeNsubjettiness=True, JetPTMin=1.0)
    jets = run_jet_test(run_generic, config, [(30.0, 0.2, 0.0, 1), (30.0, -0.2, 0.0, -1)])
    assert jets.GetEntries() == 1
    tau = [jets.At(0).Tau[i] for i in range(5)]
    assert tau[0] == pytest.approx(0.4, abs=1e-6)
    for i in range(1, 5):
        assert tau[i] == pytest.approx(0.0, abs=1e-6)

    for i in range(5):
        assert 0.0 <= tau[i] <= 1.0
    for i in range(4):
        assert tau[i] >= tau[i + 1] - 1e-6

    jets1 = run_jet_test(run_generic, config, [(30.0, 0.0, 0.0, 1)])
    for i in range(5):
        assert jets1.At(0).Tau[i] == pytest.approx(0.0, abs=1e-6)


def test_jet_output_sorted_by_pt_descending(run_generic):
    particles = [(10.0, 0.0, 0.0, 1), (40.0, 1.0, 0.0, 1), (25.0, 2.0, 0.0, 1)]
    jets = run_jet_test(run_generic, make_module_config(JetPTMin=1.0), particles)
    assert jets.GetEntries() == 3
    pts = [jets.At(i).Momentum.Pt() for i in range(3)]
    assert pts == sorted(pts, reverse=True)
    assert pts == pytest.approx([40.0, 25.0, 10.0], rel=1e-6)


def test_extreme_kinematics(run_generic):
    particles = [(1.0e4, 5.0, 0.0, 1)]
    jets = run_jet_test(run_generic, make_module_config(), particles)
    assert jets.GetEntries() == 1
    jet = jets.At(0)
    assert jet.Momentum.Pt() == pytest.approx(1.0e4, rel=1e-6)
    assert jet.Momentum.Eta() == pytest.approx(5.0, abs=1e-6)
    assert math.isfinite(jet.Momentum.E())
    assert jet.Momentum.E() == pytest.approx(1.0e4 * math.cosh(5.0), rel=1e-6)


def test_passive_area_is_cone_area(run_generic):
    config = make_module_config(AreaAlgorithm=2, ParameterR=0.5, JetPTMin=1.0)
    jets = run_jet_test(run_generic, config, [(30.0, 0.0, 0.0, 1)])
    assert jets.GetEntries() == 1
    a05 = jets.At(0).Area.E()
    assert 0.9 * math.pi * 0.25 < a05 < 1.3 * math.pi * 0.25
    config = make_module_config(AreaAlgorithm=2, ParameterR=1.0, JetPTMin=1.0)
    jets = run_jet_test(run_generic, config, [(30.0, 0.0, 0.0, 1)])
    assert jets.GetEntries() == 1
    a10 = jets.At(0).Area.E()
    assert 0.9 * math.pi * 1.0 < a10 < 1.3 * math.pi * 1.0

    assert a10 > 2.5 * a05


def test_substructure_p4_consistency(run_generic):
    config = make_module_config(
        ComputeTrimming=True, ComputePruning=True, ComputeSoftDrop=True, ParameterR=0.8, JetPTMin=1.0
    )
    jets = run_jet_test(run_generic, config, [(40.0, 0.3, 0.0, 1), (10.0, -0.3, 0.0, -1)])
    assert jets.GetEntries() == 1
    jet = jets.At(0)

    for nsub, p4 in (
        (jet.NSubJetsTrimmed, [jet.TrimmedP4[i].Pt() for i in range(5)]),
        (jet.NSubJetsPruned, [jet.PrunedP4[i].Pt() for i in range(5)]),
        (jet.NSubJetsSoftDropped, [jet.SoftDroppedP4[i].Pt() for i in range(5)]),
    ):
        populated = sum(1 for i in range(1, 5) if p4[i] > 0.0)
        assert nsub == populated

        assert p4[0] == pytest.approx(50.0, rel=1e-6)

    assert jet.SoftDroppedJet.Pt() == pytest.approx(jet.SoftDroppedP4[0].Pt(), rel=1e-9)
    assert jet.SoftDroppedSubJet1.Pt() == pytest.approx(jet.SoftDroppedP4[1].Pt(), rel=1e-9)
    assert jet.SoftDroppedSubJet2.Pt() == pytest.approx(jet.SoftDroppedP4[2].Pt(), rel=1e-9)


def test_softdrop_symmetry_cut_drops_soft_leg(run_generic):
    config = make_module_config(
        ComputeSoftDrop=True,
        SymmetryCutSoftDrop=0.1,
        BetaSoftDrop=0.0,
        R0SoftDrop=0.8,
        ParameterR=0.8,
        JetPTMin=1.0,
    )
    jets = run_jet_test(run_generic, config, [(40.0, 0.3, 0.0, 1), (3.0, -0.3, 0.0, -1)])
    assert jets.GetEntries() == 1
    jet = jets.At(0)
    assert jet.Momentum.Pt() == pytest.approx(43.0, rel=1e-6)
    assert jet.SoftDroppedJet.Pt() == pytest.approx(40.0, rel=1e-6)
    assert jet.NSubJetsSoftDropped == 0
