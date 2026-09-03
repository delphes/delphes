from itertools import combinations

import pytest
import ROOT
from .conftest import DATA_DIR, pythia8_reader_available, reader_context

data_info = [
    ("HepMC2", "test.hepmc2"),
    ("HepMC3", "test.hepmc3"),
    ("Pythia8", "test.cmnd"),
    ("STDHEP", "test.stdhep"),
]

if not pythia8_reader_available:
    data_info = [d for d in data_info if d[0] != "Pythia8"]


@pytest.mark.skipif(not pythia8_reader_available, reason="DelphesPythia8Reader not built (no Pythia8)")
def test_pythia8_reader_present_anchor():
    assert hasattr(ROOT, "DelphesPythia8Reader")


def run_reader(info, load_delphes):
    data_format, data_file = info
    with reader_context(load_delphes, data_format, DATA_DIR / data_file) as (_, factory, _, arrays, reader):
        reader.ReadEvent(factory, *arrays)
    return arrays


def extract_particles(arrays, index):
    array = arrays[index]
    size = array.GetEntries()
    return [
        {
            "PID": c.PID,
            "Status": c.Status,
            "Charge": c.Charge,
            "M1": c.M1,
            "M2": c.M2,
            "D1": c.D1,
            "D2": c.D2,
            "Px": c.Momentum.Px(),
            "Py": c.Momentum.Py(),
            "Pz": c.Momentum.Pz(),
            "E": c.Momentum.E(),
            "X": c.Position.X(),
            "Y": c.Position.Y(),
            "Z": c.Position.Z(),
            "T": c.Position.T(),
        }
        for c in (array.At(i) for i in range(size))
    ]


def sort_key(p):
    return (p["Status"], p["PID"], round(p["E"], 6))


def compare_particles(arrays_a, arrays_b, array_index):
    particles_a = sorted(extract_particles(arrays_a, array_index), key=sort_key)
    particles_b = sorted(extract_particles(arrays_b, array_index), key=sort_key)

    size_a = len(particles_a)
    size_b = len(particles_b)

    assert size_a == size_b, f"entry count {size_a} != {size_b}"

    for i, (a, b) in enumerate(zip(particles_a, particles_b)):
        assert a["PID"] == b["PID"], f"particle {i}: PID {a['PID']} != {b['PID']}"
        assert a["Status"] == b["Status"], f"particle {i}: Status {a['Status']} != {b['Status']}"
        assert a["Charge"] == b["Charge"], f"particle {i}: Charge {a['Charge']} != {b['Charge']}"
        for key in ("Px", "Py", "Pz", "E", "X", "Y", "Z", "T"):
            assert a[key] == pytest.approx(b[key], rel=1e-6), f"particle {i}: {key} {a[key]} != {b[key]}"


@pytest.mark.parametrize("a,b", list(combinations(data_info, 2)), ids=lambda p: p[0])
def test_reader_consistency(a, b, load_delphes):
    arrays_a = run_reader(a, load_delphes)
    arrays_b = run_reader(b, load_delphes)
    compare_particles(arrays_a, arrays_b, 0)
    compare_particles(arrays_a, arrays_b, 1)
    compare_particles(arrays_a, arrays_b, 2)


def in_range(idx, i1, i2):
    if i1 < 0:
        return False
    if i1 == i2 or i2 < 0:
        return idx == i1
    if i1 < i2:
        return i1 <= idx <= i2
    return idx == i1 or idx == i2


def check_mother_daughter_consistency(arrays):
    """
    Mother rules:
      M1=M2=-1 : beam particle
      M1=M2>=0 : copy of mother with changed momentum
      M1>=0, M2=-1 : decay product
      M1<M2, both>=0 : range of mothers (fragmentation if status 81-86)
      M1>M2, both>=0 : two different mothers

    Daughter rules:
      D1=D2=-1 : no daughters
      D1=D2>=0 : copy of mother with changed momentum
      D1>=0, D2=-1 : single daughter
      D1<D2, both>=0 : range of daughters
      D2<D1, both>=0 : two different daughters
    """
    particles = extract_particles(arrays, 0)
    size = len(particles)

    for i, p in enumerate(particles):
        m1, m2 = p["M1"], p["M2"]
        d1, d2 = p["D1"], p["D2"]

        for field, val in [("M1", m1), ("M2", m2), ("D1", d1), ("D2", d2)]:
            assert val == -1 or 0 <= val < size, f"particle {i}: {field}={val} out of range [-1, {size - 1}]"

        if m1 == -1 and m2 == -1:
            pass  # beam particle
        elif m1 == m2 and m1 >= 0:
            pass  # copy
        elif m1 >= 0 and m2 == -1:
            pass  # single mother
        elif m1 >= 0 and m2 >= 0 and m1 != m2:
            pass  # range or two mothers
        else:
            assert False, f"particle {i}: invalid M1={m1}, M2={m2}"

        if d1 == -1 and d2 == -1:
            pass  # no daughters
        elif d1 == d2 and d1 >= 0:
            pass  # copy
        elif d1 >= 0 and d2 == -1:
            pass  # single daughter
        elif d1 >= 0 and d2 >= 0 and d1 != d2:
            pass  # range or two daughters
        else:
            assert False, f"particle {i}: invalid D1={d1}, D2={d2}"

        # Stable particles (Status=1) should have no daughters
        if p["Status"] == 1:
            assert d1 == -1 and d2 == -1, f"stable particle {i} (Status=1) has D1={d1}, D2={d2}"

        # Incoming beams (Status=4) should have no mothers
        if p["Status"] == 4:
            assert m1 == -1 and m2 == -1, f"incoming beam {i} (Status=4) has M1={m1}, M2={m2}"

        # Verify daughter back-references
        if d1 >= 0:
            if d2 < 0 or d1 == d2:
                daughters = [d1]
            elif d1 < d2:
                daughters = range(d1, d2 + 1)
            else:
                daughters = [d1, d2]
            found = any(
                in_range(i, particles[j]["M1"], particles[j]["M2"]) for j in daughters if particles[j]["M1"] >= 0
            )
            assert found, f"particle {i} has D1={d1}, D2={d2}, but no daughter in range lists it as mother"

        # Verify mother back-references (skip beam mothers which do not list all daughters)
        if m1 >= 0:
            if m2 < 0 or m1 == m2:
                mothers = [m1]
            elif m1 < m2:
                mothers = range(m1, m2 + 1)
            else:
                mothers = [m1, m2]
            non_beam = [j for j in mothers if particles[j]["Status"] != 4]
            if non_beam:
                found = any(
                    in_range(i, particles[j]["D1"], particles[j]["D2"]) for j in non_beam if particles[j]["D1"] >= 0
                )
                assert found, f"particle {i} has M1={m1}, M2={m2}, but no mother in range lists it as daughter"


@pytest.mark.parametrize("data_info", data_info, ids=lambda p: p[0])
def test_mother_daughter_consistency(data_info, load_delphes):
    arrays = run_reader(data_info, load_delphes)
    check_mother_daughter_consistency(arrays)
