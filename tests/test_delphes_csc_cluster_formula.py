import math

import pytest
import ROOT


def test_eval_sum(load_delphes):
    formula = ROOT.DelphesCscClusterFormula()
    formula.Compile("decayR + decayZ + Ehad + Eem")
    assert formula.Eval(1.0, 2.0, 3.0, 4.0) == pytest.approx(10.0, rel=1e-9)


def test_eval_defaults(load_delphes):
    formula = ROOT.DelphesCscClusterFormula()
    formula.Compile("decayR + decayZ + Ehad + Eem")
    assert formula.Eval(1.0) == pytest.approx(1.0, rel=1e-9)


def test_eval_conditional_energy(load_delphes):
    formula = ROOT.DelphesCscClusterFormula()
    formula.Compile("(Ehad > Eem) * Ehad + (Ehad <= Eem) * Eem")
    assert formula.Eval(0.0, 0.0, 5.0, 2.0) == pytest.approx(5.0, rel=1e-9)
    assert formula.Eval(0.0, 0.0, 2.0, 5.0) == pytest.approx(5.0, rel=1e-9)


def test_whitespace_stripped(load_delphes):
    formula = ROOT.DelphesCscClusterFormula()
    formula.Compile("decayR\n + Eem")
    assert formula.Eval(1.0, 0.0, 0.0, 2.0) == pytest.approx(3.0, rel=1e-9)


SURFACE_R, SURFACE_Z, SURFACE_EHAD, SURFACE_EEM = 1.5, 2.5, 3.0, 0.5


def make_formula(expression):
    formula = ROOT.DelphesCscClusterFormula()
    formula.Compile(expression)
    return formula


@pytest.mark.parametrize(
    "expression, expected",
    [
        ("decayR + decayZ + Ehad + Eem", 7.5),
        ("Ehad / Eem", 6.0),
        ("sqrt(decayR * decayR + decayZ * decayZ)", (1.5**2 + 2.5**2) ** 0.5),
        ("abs(Eem - Ehad)", 2.5),
        ("min(Ehad, Eem)", 0.5),
        ("max(Ehad, Eem)", 3.0),
        ("exp(decayR / decayZ)", math.exp(1.5 / 2.5)),
    ],
)
def test_eval_scalar_functions(load_delphes, expression, expected):
    formula = make_formula(expression)
    assert formula.Eval(SURFACE_R, SURFACE_Z, SURFACE_EHAD, SURFACE_EEM) == pytest.approx(expected, rel=1e-9)


def test_eval_ternary_both_branches(load_delphes):
    formula = make_formula("Ehad > Eem * 2 ? 1 : 0")
    assert formula.Eval(0.0, 0.0, 3.0, 0.5) == pytest.approx(1.0, rel=1e-9)
    assert formula.Eval(0.0, 0.0, 0.5, 3.0) == pytest.approx(0.0, abs=1e-9)
