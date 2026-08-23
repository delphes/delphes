import ctypes
import struct

import pytest
import ROOT


def make_buf(n, initial=b""):
    buf = (ctypes.c_uint8 * n)()
    for i, b in enumerate(initial):
        buf[i] = b
    return buf


def as_bytes(buf, n):
    return bytes(buf[:n])


def ptr(obj):
    if isinstance(obj, ctypes.Array):
        return ctypes.cast(obj, ctypes.c_void_p)
    return ctypes.cast(ctypes.pointer(obj), ctypes.c_void_p)


def be_u32(value):
    return struct.pack(">I", value)


def new_writer(buf):
    w = ROOT.DelphesXDRWriter()
    w.SetBuffer(ptr(buf))
    return w


def new_reader(buf):
    r = ROOT.DelphesXDRReader()
    r.SetBuffer(ptr(buf))
    return r


@pytest.mark.parametrize(
    "value",
    [0x01020304, 1, -1, 0, 2**31 - 1],
)
def test_write_value_int32_big_endian_bytes(value):
    buf = make_buf(16)
    w = new_writer(buf)
    w.WriteValue(ptr(ctypes.c_int32(value)), 4)
    assert as_bytes(buf, 4) == struct.pack(">i", value)

    assert as_bytes(buf, 16)[4:] == b"\x00" * 12


@pytest.mark.parametrize(
    "value",
    [0x0102030405060708, 1, -1, 0],
)
def test_write_value_int64_big_endian_bytes(value):
    buf = make_buf(16)
    w = new_writer(buf)
    w.WriteValue(ptr(ctypes.c_int64(value)), 8)
    assert as_bytes(buf, 8) == struct.pack(">q", value)
    assert as_bytes(buf, 16)[8:] == b"\x00" * 8


@pytest.mark.parametrize(
    "value",
    [1.5, 100.0, -2.25, 0.0],
)
def test_write_value_float32_big_endian_bytes(value):
    buf = make_buf(8)
    w = new_writer(buf)
    w.WriteValue(ptr(ctypes.c_float(value)), 4)
    assert as_bytes(buf, 4) == struct.pack(">f", value)


@pytest.mark.parametrize(
    "value",
    [1.5, -1e30, 0.0],
)
def test_write_value_double_big_endian_bytes(value):
    buf = make_buf(16)
    w = new_writer(buf)
    w.WriteValue(ptr(ctypes.c_double(value)), 8)
    assert as_bytes(buf, 8) == struct.pack(">d", value)


def test_consecutive_writes_advance_offset():
    buf = make_buf(20)
    w = new_writer(buf)
    w.WriteValue(ptr(ctypes.c_int32(0x01020304)), 4)
    w.WriteValue(ptr(ctypes.c_double(1.5)), 8)
    w.WriteValue(ptr(ctypes.c_float(100.0)), 4)
    assert as_bytes(buf, 4) == struct.pack(">i", 0x01020304)
    assert as_bytes(buf, 16)[4:12] == struct.pack(">d", 1.5)
    assert as_bytes(buf, 16)[12:16] == struct.pack(">f", 100.0)


def test_writer_set_offset_starts_at_offset():
    buf = make_buf(16, initial=b"\xaa" * 16)
    w = new_writer(buf)
    w.SetOffset(8)
    w.WriteValue(ptr(ctypes.c_int32(0x01020304)), 4)
    assert as_bytes(buf, 8) == b"\xaa" * 8
    assert as_bytes(buf, 12)[8:12] == struct.pack(">i", 0x01020304)
    assert as_bytes(buf, 16)[12:] == b"\xaa" * 4


def test_writer_set_buffer_resets_offset_to_zero():
    buf = make_buf(8)
    w = new_writer(buf)
    w.WriteValue(ptr(ctypes.c_int32(0x11111111)), 4)
    w.SetOffset(4)
    w.SetBuffer(ptr(buf))
    w.WriteValue(ptr(ctypes.c_int32(0x22222222)), 4)
    assert as_bytes(buf, 4) == struct.pack(">i", 0x22222222)


def test_write_raw_is_noop_in_buffer_mode():
    buf = make_buf(8, initial=b"\xde" * 8)
    w = new_writer(buf)
    src = (ctypes.c_uint8 * 4)(1, 2, 3, 4)
    w.WriteRaw(ptr(src), 4)
    assert as_bytes(buf, 8) == b"\xde" * 8


@pytest.mark.parametrize(
    "wire, value",
    [
        (b"\x01\x02\x03\x04", 0x01020304),
        (b"\x00\x00\x00\x01", 1),
        (b"\xff\xff\xff\xff", -1),
        (b"\x7f\xff\xff\xff", 2**31 - 1),
    ],
)
def test_read_value_int32_hand_crafted(wire, value):
    buf = make_buf(8, initial=wire)
    r = new_reader(buf)
    out = ctypes.c_int32(0)
    r.ReadValue(ptr(out), 4)
    assert out.value == value


@pytest.mark.parametrize(
    "wire, value",
    [
        (b"\x01\x02\x03\x04\x05\x06\x07\x08", 0x0102030405060708),
        (b"\x00" * 7 + b"\x01", 1),
        (b"\xff" * 8, -1),
    ],
)
def test_read_value_int64_hand_crafted(wire, value):
    buf = make_buf(12, initial=wire)
    r = new_reader(buf)
    out = ctypes.c_int64(0)
    r.ReadValue(ptr(out), 8)
    assert out.value == value


@pytest.mark.parametrize(
    "wire, value",
    [
        (b"\x3f\xc0\x00\x00", 1.5),
        (b"\x42\xc8\x00\x00", 100.0),
        (b"\xc0\x10\x00\x00", -2.25),
    ],
)
def test_read_value_float32_hand_crafted(wire, value):
    buf = make_buf(8, initial=wire)
    r = new_reader(buf)
    out = ctypes.c_float(0.0)
    r.ReadValue(ptr(out), 4)
    assert out.value == value


@pytest.mark.parametrize(
    "wire, value",
    [
        (b"\x3f\xf8" + b"\x00" * 6, 1.5),
        (b"\x00" * 8, 0.0),
        (struct.pack(">d", 0.1), 0.1),
    ],
)
def test_read_value_double_hand_crafted(wire, value):
    buf = make_buf(12, initial=wire)
    r = new_reader(buf)
    out = ctypes.c_double(0.0)
    r.ReadValue(ptr(out), 8)
    assert out.value == value


def test_reader_set_offset_reads_correct_slice():
    prefix = b"\x12\x34\x56\x78"
    buf = make_buf(8, initial=prefix + be_u32(0x01020304))
    r = new_reader(buf)
    r.SetOffset(4)
    out = ctypes.c_int32(0)
    r.ReadValue(ptr(out), 4)
    assert out.value == 0x01020304


def test_reader_set_buffer_resets_offset_to_zero():
    buf = make_buf(8, initial=be_u32(0x11111111) + be_u32(0x22222222))
    r = new_reader(buf)
    r.ReadValue(ptr(ctypes.c_int32()), 4)
    r.SetBuffer(ptr(buf))
    out = ctypes.c_int32(0)
    r.ReadValue(ptr(out), 4)
    assert out.value == 0x11111111


def test_consecutive_reads_advance_offset():
    buf = make_buf(16, initial=be_u32(0x01020304) + struct.pack(">d", 1.5) + be_u32(0xAABBCCDD))
    r = new_reader(buf)
    out32 = ctypes.c_int32(0)
    r.ReadValue(ptr(out32), 4)
    out64 = ctypes.c_double(0.0)
    r.ReadValue(ptr(out64), 8)
    assert out32.value == 0x01020304
    assert out64.value == 1.5


def test_read_raw_is_noop_in_buffer_mode():
    buf = make_buf(8, initial=be_u32(0x01020304))
    r = new_reader(buf)
    out = (ctypes.c_uint8 * 4)()
    r.ReadRaw(ptr(out), 4)
    assert bytes(out) == b"\x00" * 4


def test_read_string_basic_no_trailing_nul():
    buf = make_buf(12, initial=be_u32(5) + b"hello")
    r = new_reader(buf)
    dst = (ctypes.c_uint8 * 8)(*([0xEE] * 8))
    r.ReadString(ptr(dst), 8)
    assert bytes(dst[:5]) == b"hello"

    assert dst[5] == 0xEE
    assert dst[7] == 0xEE


def test_read_string_stored_trailing_nul():
    buf = make_buf(12, initial=be_u32(6) + b"hello\x00")
    r = new_reader(buf)
    dst = (ctypes.c_uint8 * 8)(*([0xEE] * 8))
    r.ReadString(ptr(dst), 8)
    assert bytes(dst[:6]) == b"hello\x00"
    assert dst[6] == 0xEE


def test_read_string_embedded_nul_not_truncated():
    buf = make_buf(12, initial=be_u32(4) + b"a\x00bc")
    r = new_reader(buf)
    dst = (ctypes.c_uint8 * 8)(*([0xEE] * 8))
    r.ReadString(ptr(dst), 8)
    assert bytes(dst[:4]) == b"a\x00bc"


def test_read_string_empty():
    buf = make_buf(8, initial=be_u32(0) + be_u32(0x01020304))
    r = new_reader(buf)
    dst = (ctypes.c_uint8 * 4)(*([0xEE] * 4))
    r.ReadString(ptr(dst), 8)
    assert bytes(dst) == b"\xee" * 4

    out = ctypes.c_int32(0)
    r.ReadValue(ptr(out), 4)
    assert out.value == 0x01020304


def test_read_string_max_size_shorter_clamps_copy_and_offset():
    buf = make_buf(16, initial=be_u32(8) + b"abcdefgh" + be_u32(0x0000CAFE))
    r = new_reader(buf)
    dst = (ctypes.c_uint8 * 8)(*([0xEE] * 8))
    r.ReadString(ptr(dst), 3)
    assert bytes(dst[:3]) == b"abc"
    assert bytes(dst[3:]) == b"\xee" * 5

    out = ctypes.c_int32(0)
    r.ReadValue(ptr(out), 4)
    assert out.value == 0x64656667


def test_read_string_max_size_longer_than_stored():
    buf = make_buf(16, initial=be_u32(5) + b"hello" + be_u32(0x01020304))
    r = new_reader(buf)
    dst = (ctypes.c_uint8 * 12)(*([0xEE] * 12))
    r.ReadString(ptr(dst), 12)
    assert bytes(dst[:5]) == b"hello"
    assert bytes(dst[5:]) == b"\xee" * 7
    out = ctypes.c_int32(0)
    r.ReadValue(ptr(out), 4)
    assert out.value == 0x01020304


def test_mixed_sequence_roundtrip():
    n_particles = 3
    px, py, pz = -12.5, 3.25, 1000.0
    tag = 0x01020304
    string_bytes = b"pile-up event"

    buf = make_buf(64)
    w = new_writer(buf)
    w.WriteValue(ptr(ctypes.c_int32(n_particles)), 4)
    w.WriteValue(ptr(ctypes.c_float(px)), 4)
    w.WriteValue(ptr(ctypes.c_float(py)), 4)
    w.WriteValue(ptr(ctypes.c_double(pz)), 8)
    w.WriteValue(ptr(ctypes.c_uint32(tag)), 4)

    head = 24
    for i, b in enumerate(be_u32(len(string_bytes)) + string_bytes):
        buf[head + i] = b

    r = new_reader(buf)
    count = ctypes.c_int32(0)
    r.ReadValue(ptr(count), 4)
    x = ctypes.c_float(0.0)
    r.ReadValue(ptr(x), 4)
    y = ctypes.c_float(0.0)
    r.ReadValue(ptr(y), 4)
    z = ctypes.c_double(0.0)
    r.ReadValue(ptr(z), 8)
    out_tag = ctypes.c_uint32(0)
    r.ReadValue(ptr(out_tag), 4)
    dst = (ctypes.c_uint8 * 32)()
    r.ReadString(ptr(dst), 32)

    assert count.value == n_particles
    assert x.value == px
    assert y.value == py
    assert z.value == pz
    assert out_tag.value == tag
    assert bytes(dst[: len(string_bytes)]) == string_bytes
