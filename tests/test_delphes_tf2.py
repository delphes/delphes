import ctypes
import math

import pytest
import ROOT


def make_tf2(expression):
    func = ROOT.DelphesTF2()
    func.Compile(expression)
    return func


def test_eval(load_delphes):
    func = make_tf2("z*z + t*t")
    assert func.Eval(3.0, 4.0) == pytest.approx(25.0, rel=1e-9)


def test_eval_linear(load_delphes):
    func = make_tf2("2*z + 3*t + 1")
    assert func.Eval(1.0, 2.0) == pytest.approx(9.0, rel=1e-9)


def test_get_random2_within_range(load_delphes):
    func = make_tf2("z*z + t*t")
    func.SetRange(-2.0, -2.0, 2.0, 2.0)
    ROOT.gRandom.SetSeed(42)
    x = ctypes.c_double()
    y = ctypes.c_double()
    for _ in range(200):
        func.GetRandom2(x, y)
        assert -2.0 <= x.value <= 2.0
        assert -2.0 <= y.value <= 2.0


def test_get_random2_deterministic(load_delphes):
    func = make_tf2("z*z + t*t")
    func.SetRange(-2.0, -2.0, 2.0, 2.0)
    ROOT.gRandom.SetSeed(42)
    x1 = ctypes.c_double()
    y1 = ctypes.c_double()
    func.GetRandom2(x1, y1)
    ROOT.gRandom.SetSeed(42)
    x2 = ctypes.c_double()
    y2 = ctypes.c_double()
    func.GetRandom2(x2, y2)
    assert x1.value == x2.value
    assert y1.value == y2.value


def sample_random2(func, n, seed=42):
    ROOT.gRandom.SetSeed(seed)
    xs, ys = [], []
    for _ in range(n):
        x = ctypes.c_double()
        y = ctypes.c_double()
        func.GetRandom2(x, y)
        xs.append(x.value)
        ys.append(y.value)
    return xs, ys


def test_eval_nontrivial_2d_values(load_delphes):
    assert make_tf2("z * t").Eval(3.0, 4.0) == pytest.approx(12.0, rel=1e-9)
    assert make_tf2("z^2 - t^2").Eval(3.0, 4.0) == pytest.approx(-7.0, rel=1e-9)
    assert make_tf2("sin(z) * cos(t)").Eval(3.0, 4.0) == pytest.approx(math.sin(3.0) * math.cos(4.0), rel=1e-9)


def test_get_random2_non_square_range(load_delphes):
    func = make_tf2("z*z + t*t")
    func.SetRange(0.0, -1.0, 10.0, 1.0)
    xs, ys = sample_random2(func, 300)
    assert all(0.0 <= x <= 10.0 for x in xs)
    assert all(-1.0 <= y <= 1.0 for y in ys)
    assert min(xs) < max(xs)
    assert min(ys) < max(ys)


def test_get_random2_two_argument_range_only_x_axis(load_delphes):
    func = make_tf2("z*z + t*t")
    func.SetRange(0.0, -1.0, 10.0, 1.0)
    func.SetRange(2.0, 5.0)
    xs, ys = sample_random2(func, 300)
    assert all(2.0 <= x <= 5.0 for x in xs)
    assert all(-1.0 <= y <= 1.0 for y in ys)


def test_get_random_is_stub(load_delphes):
    func = make_tf2("z*z + t*t")
    func.SetRange(0.0, 0.0, 10.0, 10.0)
    ROOT.gRandom.SetSeed(1)
    for _ in range(10):
        assert func.GetRandom() == 0.0
