import math

import pytest
import ROOT


def make_formula(expression):
    formula = ROOT.DelphesFormula()
    formula.Compile(expression)
    return formula


SURFACE_PT, SURFACE_ETA, SURFACE_PHI, SURFACE_ENERGY = 2.0, 0.5, 0.25, 16.0


@pytest.mark.parametrize(
    "expression, expected",
    [
        ("log(pt)", math.log(2.0)),
        ("log10(pt)", math.log10(2.0)),
        ("log1p(pt)", math.log1p(2.0)),
        ("log2(pt)", 1.0),
        ("exp(pt)", math.exp(2.0)),
        ("sqrt(pt)", math.sqrt(2.0)),
        ("sin(pt)", math.sin(2.0)),
        ("cos(pt)", math.cos(2.0)),
        ("tan(pt)", math.tan(2.0)),
        ("asin(eta)", math.asin(0.5)),
        ("acos(eta)", math.acos(0.5)),
        ("atan(phi)", math.atan(0.25)),
        ("atan2(energy, pt)", math.atan2(16.0, 2.0)),
        ("sinh(pt)", math.sinh(2.0)),
        ("cosh(pt)", math.cosh(2.0)),
        ("tanh(pt)", math.tanh(2.0)),
        ("min(energy, pt)", 2.0),
        ("max(energy, pt)", 16.0),
        ("pow(energy, 0.5)", 4.0),
        ("pt ** 2", 4.0),
        ("floor(energy / 5)", 3.0),
        ("trunc(energy / 5)", 3.0),
        ("abs(1 - pt)", 1.0),
    ],
)
def test_eval_scalar_functions(load_delphes, expression, expected):
    formula = make_formula(expression)
    assert formula.Eval(SURFACE_PT, SURFACE_ETA, SURFACE_PHI, SURFACE_ENERGY) == pytest.approx(expected, rel=1e-9)


def test_eval_operator_precedence(load_delphes):
    assert make_formula("1 + 2 * pt").Eval(2.0) == pytest.approx(5.0, rel=1e-9)
    assert make_formula("energy / pt * pt").Eval(2.0, 0.0, 0.0, 16.0) == pytest.approx(16.0, rel=1e-9)
    assert make_formula("(1 + 2) * pt").Eval(2.0) == pytest.approx(6.0, rel=1e-9)

    assert make_formula("-pt + 2").Eval(2.0) == pytest.approx(0.0, abs=1e-9)


def test_eval_ternary_both_branches(load_delphes):
    formula = make_formula("pt < 5 ? 10 : 20")
    assert formula.Eval(4.0) == pytest.approx(10.0, rel=1e-9)
    assert formula.Eval(6.0) == pytest.approx(20.0, rel=1e-9)


def test_eval_ctg_theta_candidate_parameter(load_delphes):
    _, factory, _ = load_delphes({})
    formula = make_formula("d0 + dz + ctgTheta + radius + density")
    candidate = factory.NewCandidate()
    candidate.D0 = 1.0
    candidate.DZ = 2.0
    candidate.CtgTheta = 3.0
    candidate.Position.SetPtEtaPhiE(4.0, 0.0, 0.0, 4.0)
    candidate.ParticleDensity = 5.0
    assert formula.Eval(0.0, 0.0, 0.0, 0.0, candidate) == pytest.approx(15.0, rel=1e-9)

    assert formula.Eval(0.0, 0.0, 0.0, 0.0) == pytest.approx(0.0, abs=1e-9)


def test_eval_division_by_zero_is_infinity(load_delphes):
    value = make_formula("energy / pt").Eval(0.0, 0.0, 0.0, 16.0)
    assert math.isinf(value)
    assert value > 0.0


def test_eval_domain_errors_are_nan_or_infinity(load_delphes):
    assert math.isnan(make_formula("sqrt(-pt)").Eval(2.0))
    assert make_formula("log(0.0 * pt)").Eval(2.0) == pytest.approx(-math.inf)
    assert math.isnan(make_formula("asin(energy)").Eval(2.0, 0.0, 0.0, 16.0))


def test_eval_nested_expression_precision(load_delphes):
    formula = make_formula("(pt + eta) * (phi - energy) / (pt * pt + 1) + sqrt(energy) * 0.5")
    assert formula.Eval(2.0, 0.5, 0.25, 16.0) == pytest.approx(-5.875, rel=1e-12)


def test_eval_basic_sum(load_delphes):
    formula = make_formula("pt + eta + phi + energy")
    assert formula.Eval(1.0, 2.0, 3.0, 4.0) == pytest.approx(10.0, rel=1e-9)


def test_eval_pt_only_with_defaults(load_delphes):
    formula = make_formula("pt * pt")
    assert formula.Eval(3.0) == pytest.approx(9.0, rel=1e-9)

    assert formula.Eval(3.0, 0.5, 0.0, 100.0) == pytest.approx(9.0, rel=1e-9)


def test_eval_conditional_on_pt(load_delphes):
    formula = make_formula("(pt > 10) * 2 + (pt <= 10) * 0.5")
    assert formula.Eval(20.0) == pytest.approx(2.0, rel=1e-9)
    assert formula.Eval(10.0) == pytest.approx(0.5, rel=1e-9)
    assert formula.Eval(5.0) == pytest.approx(0.5, rel=1e-9)


def test_eval_abs(load_delphes):
    formula = make_formula("abs(pt - 10)")
    assert formula.Eval(7.0) == pytest.approx(3.0, rel=1e-9)
    assert formula.Eval(13.0) == pytest.approx(3.0, rel=1e-9)


def test_eval_division(load_delphes):
    formula = make_formula("energy / pt")
    assert formula.Eval(2.0, 0.0, 0.0, 6.0) == pytest.approx(3.0, rel=1e-9)


def test_eval_with_candidate_parameters(load_delphes):
    _, factory, _ = load_delphes({})
    formula = make_formula("d0 + dz + radius + density")
    candidate = factory.NewCandidate()
    candidate.D0 = 1.0
    candidate.DZ = 2.0
    candidate.Position.SetPtEtaPhiE(4.0, 0.0, 0.0, 4.0)
    candidate.ParticleDensity = 5.0
    assert formula.Eval(0.0, 0.0, 0.0, 0.0, candidate) == pytest.approx(12.0, rel=1e-9)

    assert formula.Eval(0.0, 0.0, 0.0, 0.0) == pytest.approx(0.0, abs=1e-9)


def test_whitespace_and_backslash_stripped(load_delphes):
    formula = make_formula("pt\n + 1")
    assert formula.Eval(2.0) == pytest.approx(3.0, rel=1e-9)
    formula = make_formula("pt \\\t + 2")
    assert formula.Eval(1.0) == pytest.approx(3.0, rel=1e-9)
