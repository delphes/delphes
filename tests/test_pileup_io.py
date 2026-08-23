import ctypes
import gc
import math
import struct

import pytest
import ROOT
from conftest import MINBIAS_FILE

EVENTS = [
    [(211, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0), (-13, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0)],
    [(22, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0)],
]


def write_pileup_events(path, events):
    writer = ROOT.DelphesPileUpWriter(str(path))
    for event in events:
        for pid, x, y, z, t, px, py, pz, e in event:
            writer.WriteParticle(pid, x, y, z, t, px, py, pz, e)
        writer.WriteEntry()
    writer.WriteIndex()

    del writer
    gc.collect()


def fresh_read_fields():
    pid = ctypes.c_int32()
    fields = [ctypes.c_float() for _ in range(8)]
    return pid, fields


def read_event(reader, entry, expected):
    assert reader.ReadEntry(entry)
    pid, fields = fresh_read_fields()
    for pid_expected, *values in expected:
        assert reader.ReadParticle(pid, *fields)
        assert pid.value == pid_expected
        for field, value in zip(fields, values):
            assert field.value == pytest.approx(value, rel=1e-6)

    assert not reader.ReadParticle(pid, *fields)


def test_write_read_roundtrip(tmp_path):
    path = tmp_path / "roundtrip.pileup"
    write_pileup_events(path, EVENTS)
    reader = ROOT.DelphesPileUpReader(str(path))
    assert reader.GetEntries() == 2
    read_event(reader, 0, EVENTS[0])
    read_event(reader, 1, EVENTS[1])

    read_event(reader, 0, EVENTS[0])


def test_read_out_of_range_entry(tmp_path):
    path = tmp_path / "roundtrip.pileup"
    write_pileup_events(path, EVENTS)
    reader = ROOT.DelphesPileUpReader(str(path))
    assert not reader.ReadEntry(2)
    assert not reader.ReadEntry(100)


def test_real_minbias_pileup_parses():
    reader = ROOT.DelphesPileUpReader(MINBIAS_FILE)
    entries = reader.GetEntries()
    assert entries > 0
    assert reader.ReadEntry(0)
    pid, fields = fresh_read_fields()
    count = 0
    while reader.ReadParticle(pid, *fields):
        count += 1
        px, py, _, e = (f.value for f in fields[4:])
        assert e > 0.0
        assert px * px + py * py > 0.0
        assert pid.value != 0
    assert count > 0


def float32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def read_particle_values(reader):
    pid, fields = fresh_read_fields()
    assert reader.ReadParticle(pid, *fields)
    return [pid.value] + [f.value for f in fields]


def test_single_particle_event(tmp_path):
    path = tmp_path / "single.pileup"
    write_pileup_events(path, [[(211, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)]])

    reader = ROOT.DelphesPileUpReader(str(path))
    assert reader.GetEntries() == 1
    assert reader.ReadEntry(0)
    assert read_particle_values(reader) == [211, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

    pid, fields = fresh_read_fields()
    assert not reader.ReadParticle(pid, *fields)


def test_empty_event_roundtrip(tmp_path):
    path = tmp_path / "empty.pileup"
    write_pileup_events(
        path,
        [
            [(22, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0), (13, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)],
            [],
            [(11, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)],
        ],
    )

    reader = ROOT.DelphesPileUpReader(str(path))
    assert reader.GetEntries() == 3
    assert reader.ReadEntry(0)
    assert read_particle_values(reader)[0] == 22
    assert read_particle_values(reader)[0] == 13

    assert reader.ReadEntry(1)
    pid, fields = fresh_read_fields()
    assert not reader.ReadParticle(pid, *fields)

    assert reader.ReadEntry(2)
    assert read_particle_values(reader) == [11, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]


def test_float32_precision_boundary(tmp_path):
    path = tmp_path / "float32.pileup"
    exact = 0.1
    large = 1.0e30
    assert float32(exact) != exact
    assert float32(large) != large

    write_pileup_events(path, [[(211, exact, 2.0, large, 4.0, 5.0, 6.0, 7.0, 8.0)]])

    reader = ROOT.DelphesPileUpReader(str(path))
    assert reader.ReadEntry(0)
    values = read_particle_values(reader)
    assert values[1] == float32(exact)
    assert values[1] != exact
    assert values[3] == float32(large)
    assert values[3] != large

    assert values[2] == 2.0


def test_negative_and_neg_zero_signs(tmp_path):
    path = tmp_path / "signs.pileup"
    write_pileup_events(path, [[(-211, -1.5, -0.0, -3.25, -4.0, -5.0, -6.0, -7.0, -8.0)]])

    reader = ROOT.DelphesPileUpReader(str(path))
    assert reader.ReadEntry(0)
    values = read_particle_values(reader)
    assert values[0] == -211
    assert values[1] == -1.5
    assert values[2] == -0.0
    assert math.copysign(1.0, values[2]) == -1.0
    assert values[3] == -3.25
    assert values[5] == -5.0
    assert values[8] == -8.0


def test_many_entries_out_of_order_seek(tmp_path):
    n_per_entry = 2000
    path = tmp_path / "many.pileup"
    write_pileup_events(
        path,
        [
            [(211 + i, entry * 1000.0 + i, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0) for i in range(n_per_entry)]
            for entry in range(3)
        ],
    )

    reader = ROOT.DelphesPileUpReader(str(path))
    assert reader.GetEntries() == 3
    for entry in (2, 0, 1):
        assert reader.ReadEntry(entry)
        first = read_particle_values(reader)
        assert first[0] == 211
        assert first[1] == entry * 1000.0
        for _ in range(n_per_entry - 2):
            read_particle_values(reader)
        last = read_particle_values(reader)
        assert last[0] == 211 + n_per_entry - 1
        assert last[1] == entry * 1000.0 + n_per_entry - 1
        pid, fields = fresh_read_fields()
        assert not reader.ReadParticle(pid, *fields)
