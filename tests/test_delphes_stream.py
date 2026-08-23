import ctypes
import math

import pytest
import ROOT


def make_stream(buffer):
    return ROOT.DelphesStream(buffer)


def test_read_double_and_int(load_delphes):
    stream = ROOT.DelphesStream(b"1.5 2 abc")
    value = ctypes.c_double()
    assert stream.ReadDbl(value)
    assert value.value == 1.5
    integer = ctypes.c_int()
    assert stream.ReadInt(integer)
    assert integer.value == 2
    assert stream.FindStr(b"abc")


def test_find_chr(load_delphes):
    stream = ROOT.DelphesStream(b"1.5;2.5")
    value = ctypes.c_double()
    assert stream.ReadDbl(value)
    assert value.value == 1.5
    assert stream.FindChr(ord(";"))
    assert stream.ReadDbl(value)
    assert value.value == 2.5


def test_read_double_malformed(load_delphes):
    stream = ROOT.DelphesStream(b"abc")
    value = ctypes.c_double()
    assert not stream.ReadDbl(value)


def test_read_int_malformed(load_delphes):
    stream = ROOT.DelphesStream(b"xyz")
    integer = ctypes.c_int()
    assert not stream.ReadInt(integer)


def test_find_str_missing(load_delphes):
    stream = ROOT.DelphesStream(b"1 2 3")
    assert not stream.FindStr(b"nonexistent")


def test_leading_whitespace_and_signs(load_delphes):
    stream = ROOT.DelphesStream(b"  -3.5e2 +7 0.25")
    value = ctypes.c_double()
    assert stream.ReadDbl(value)
    assert value.value == -350.0
    integer = ctypes.c_int()
    assert stream.ReadInt(integer)
    assert integer.value == 7
    assert stream.ReadDbl(value)
    assert value.value == 0.25


@pytest.mark.parametrize(
    "text, expected",
    [
        ("2147483647", 2147483647),
        ("-2147483648", -2147483648),
    ],
)
def test_read_int_range_extremes(load_delphes, text, expected):
    integer = ctypes.c_int()
    assert make_stream(text.encode()).ReadInt(integer)
    assert integer.value == expected


def test_read_int_outrange_clamped_first_in_process(load_delphes):
    integer = ctypes.c_int()
    assert make_stream(b"2147483648").ReadInt(integer)
    assert integer.value == 2147483647
    integer = ctypes.c_int()
    assert make_stream(b"-2147483649").ReadInt(integer)
    assert integer.value == -2147483648


@pytest.mark.parametrize(
    "text, expected",
    [
        (".5", 0.5),
        ("5.", 5.0),
        ("1e-3", 0.001),
        ("1e+003", 1000.0),
        ("+0.0", 0.0),
        ("1e400", math.inf),
        ("1e-400", 0.0),
    ],
)
def test_read_double_edge_numbers(load_delphes, text, expected):
    value = ctypes.c_double()
    assert make_stream(text.encode()).ReadDbl(value)
    assert value.value == expected


def test_read_double_negative_zero(load_delphes):
    value = ctypes.c_double()
    assert make_stream(b"-0.0").ReadDbl(value)
    assert value.value == 0.0
    assert math.copysign(1.0, value.value) == -1.0


def test_read_at_end_returns_false_without_move(load_delphes):
    stream = make_stream(b"1.5")
    value = ctypes.c_double()
    assert stream.ReadDbl(value)
    assert value.value == 1.5
    sentinel = ctypes.c_double(99.0)
    assert not stream.ReadDbl(sentinel)
    assert sentinel.value == 0.0
    integer = ctypes.c_int(7)
    assert not stream.ReadInt(integer)
    assert integer.value == 0


def test_read_malformed_does_not_move_buffer(load_delphes):
    stream = make_stream(b"abc15")
    value = ctypes.c_double(99.0)
    assert not stream.ReadDbl(value)
    assert value.value == 0.0
    assert stream.FindChr(ord("1"))
    assert stream.ReadDbl(value)
    assert value.value == 5.0


def test_empty_buffer(load_delphes):
    stream = make_stream(b"")
    value = ctypes.c_double()
    integer = ctypes.c_int()
    assert not stream.ReadDbl(value)
    assert not stream.ReadInt(integer)
    assert not stream.FindStr(b"a")
    assert not stream.FindChr(ord("a"))


def test_whitespace_only_buffer(load_delphes):
    value = ctypes.c_double()
    assert not make_stream(b"   ").ReadDbl(value)


def test_find_chr_first_and_last_character(load_delphes):
    assert make_stream(b"abc").FindChr(ord("a"))

    assert make_stream(b"1.5").FindChr(ord("5"))


def test_find_chr_zero_finds_terminator(load_delphes):
    assert make_stream(b"abc").FindChr(0)


def test_find_str_first_and_last(load_delphes):
    stream = make_stream(b"xyzabc")
    assert stream.FindStr(b"xyz")
    assert stream.FindStr(b"abc")
    stream = make_stream(b"abcxyz")
    assert stream.FindStr(b"xyz")
