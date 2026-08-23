import math

import pytest
from conftest import make_candidate, make_config


def test_scale_by_factor(run_module):
    config = make_config("EnergyScale", ScaleFormula=0.5)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 25.0


def test_scale_by_one(run_module):
    config = make_config("EnergyScale", ScaleFormula=1.0)
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}])
    assert output.GetEntries() == 1
    assert output.At(0).Momentum.Pt() == 50.0


def test_scale_formula_by_pt(run_module):
    config = make_config("EnergyScale", ScaleFormula="{(pt > 10) * 0.8 + (pt <= 10) * 0.5}")
    output = run_module(config, [{"pt": 50.0, "eta": 0.5}, {"pt": 5.0, "eta": 0.5}])
    assert output.GetEntries() == 2
    assert output.At(0).Momentum.Pt() == 40.0
    assert output.At(1).Momentum.Pt() == 2.5


def test_preserves_eta_phi(run_module):
    config = make_config("EnergyScale", ScaleFormula=0.5)
    output = run_module(config, [{"pt": 50.0, "eta": 1.2, "phi": 0.7}])
    smeared = output.At(0)
    assert smeared.Momentum.Eta() == 1.2
    assert smeared.Momentum.Phi() == pytest.approx(0.7, abs=1e-6)


def test_empty_input(run_module):
    output = run_module(make_config("EnergyScale", ScaleFormula=0.5), [])
    assert output.GetEntries() == 0


def test_massive_particle_full_p4_linearity(run_generic):
    pt, eta, phi, m = 50.0, 0.5, 0.3, 10.0
    scale = 2.0

    def setup(module, factory):
        arr = module.ExportArray("inputParticles")
        c = make_candidate(factory, pt, eta, phi, pid=211, charge=1)
        c.Momentum.SetPtEtaPhiM(pt, eta, phi, m)
        arr.Add(c)

    config = make_config("EnergyScale", ScaleFormula=scale)
    output = run_generic(config, setup=setup)
    assert output.GetEntries() == 1
    mom = output.At(0).Momentum
    e_in = math.sqrt(pt * pt * math.cosh(eta) ** 2 + m * m)
    assert mom.Pt() == pytest.approx(scale * pt, rel=1e-9)
    assert mom.Py() == pytest.approx(scale * pt * math.sin(phi), rel=1e-9)
    assert mom.Pz() == pytest.approx(scale * pt * math.sinh(eta), rel=1e-9)
    assert mom.E() == pytest.approx(scale * e_in, rel=1e-9)

    assert mom.M() == pytest.approx(scale * m, rel=1e-9)

    assert mom.M2() == pytest.approx(mom.E() ** 2 - mom.Px() ** 2 - mom.Py() ** 2 - mom.Pz() ** 2, rel=1e-9)
