import math

import pytest
import ROOT


def test_eval_sum(load_delphes):
    formula = ROOT.DelphesCylindricalFormula()
    formula.Compile("rho + phi + z")
    assert formula.Eval(1.0, 2.0, 3.0) == pytest.approx(6.0, rel=1e-9)


def test_eval_defaults(load_delphes):
    formula = ROOT.DelphesCylindricalFormula()
    formula.Compile("rho + phi + z")

    assert formula.Eval(1.0) == pytest.approx(1.0, rel=1e-9)


def test_eval_pythagoras(load_delphes):
    formula = ROOT.DelphesCylindricalFormula()
    formula.Compile("rho*rho + z*z")
    assert formula.Eval(3.0, 0.0, 4.0) == pytest.approx(25.0, rel=1e-9)


def test_whitespace_stripped(load_delphes):
    formula = ROOT.DelphesCylindricalFormula()
    formula.Compile("rho\n + z")
    assert formula.Eval(1.0, 0.0, 2.0) == pytest.approx(3.0, rel=1e-9)


SURFACE_R, SURFACE_PHI, SURFACE_Z = 3.0, 0.1, 4.0


def make_formula(expression):
    formula = ROOT.DelphesCylindricalFormula()
    formula.Compile(expression)
    return formula


@pytest.mark.parametrize(
    "expression, expected",
    [
        ("rho * rho + z * z", 25.0),
        ("sqrt(rho * rho) * z", 12.0),
        ("floor(rho + 0.5)", 3.0),
        ("exp(z / rho)", math.exp(4.0 / 3.0)),
        ("pow(rho, 2) + sin(phi) * z", 9.0 + math.sin(0.1) * 4.0),
        ("abs(z - rho)", 1.0),
        ("min(rho, z)", 3.0),
        ("max(rho, z)", 4.0),
        ("rho * phi + z", 4.3),
    ],
)
def test_eval_scalar_functions(load_delphes, expression, expected):
    formula = make_formula(expression)
    assert formula.Eval(SURFACE_R, SURFACE_PHI, SURFACE_Z) == pytest.approx(expected, rel=1e-9)


def test_eval_conditional_both_branches(load_delphes):
    formula = make_formula("rho < 5 ? z : 0")
    assert formula.Eval(4.0, 0.0, 7.0) == pytest.approx(7.0, rel=1e-9)
    assert formula.Eval(6.0, 0.0, 7.0) == pytest.approx(0.0, abs=1e-9)


def test_eval_boundaries_zero_r_phi_z(load_delphes):
    formula = make_formula("rho * rho + z * z")

    assert formula.Eval(0.0, 1.0, 5.0) == pytest.approx(25.0, rel=1e-9)
    assert formula.Eval(3.0, 1.0, 0.0) == pytest.approx(9.0, rel=1e-9)
    assert formula.Eval(3.0, 0.1, -4.0) == pytest.approx(25.0, rel=1e-9)
