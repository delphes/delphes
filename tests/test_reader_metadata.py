import contextlib
import math

import pytest
import ROOT
from .conftest import DATA_DIR, check_event_fields, read_tree_branch, reader_context

HEPMCE2_EVENT_1 = {
    "number": 0,
    "mpi": 2,
    "scale": 1.8177109999999999e2,
    "alpha_qcd": 1.1661000000000001e-1,
    "alpha_qed": 7.5467709999999999e-3,
    "process_id": 9999,
    "weight": 5.0503900000000004e2,
    "cross_section": 5.0503900000000004e2,
    "cross_section_error": 5.0503900000000004e2,
    "id1": 21,
    "id2": 21,
    "x1": 3.6198262818461535e-1,
    "x2": 2.1623346607692307e-3,
    "scale_pdf": 1.8177109999999999e2,
    "pdf1": 0.0,
    "pdf2": 0.0,
}


HEPMCE3_EVENT_1 = {
    "number": 0,
    "mpi": 2,
    "scale": 181.7711,
    "alpha_qcd": 0.11661,
    "alpha_qed": 0.007546771,
    "process_id": 9999,
    "weight": 5.0503900000000004411049e2,
    "cross_section": 5.05039000e2,
    "cross_section_error": 5.05039000e2,
    "id1": 21,
    "id2": 21,
    "x1": 3.61982628e-1,
    "x2": 2.16233466e-3,
    "scale_pdf": 1.81771100e2,
    "pdf1": 0.0,
    "pdf2": 0.0,
}


@contextlib.contextmanager
def make_reader(format, data_file, load_delphes):
    with reader_context(load_delphes, format, DATA_DIR / data_file) as (_, factory, _, arrays, reader):
        assert reader.ReadEvent(factory, *arrays)
        yield reader


def test_hepmc2_event_metadata(load_delphes, tmp_path):
    with make_reader("HepMC2", "test.hepmc2", load_delphes) as reader:

        def build(module, writer):
            branch = writer.NewBranch("Event", ROOT.HepMCEvent.Class())
            reader.AnalyzeEvent(branch, 0)

        tree_reader, fin = read_tree_branch(load_delphes, tmp_path, build)
        branch = tree_reader.UseBranch("Event")
        tree_reader.ReadEntry(0)
        assert branch.GetEntries() == 1
        check_event_fields(branch.At(0), HEPMCE2_EVENT_1)
        fin.Close()


def test_hepmc2_analyze_weight(load_delphes, tmp_path):
    with make_reader("HepMC2", "test.hepmc2", load_delphes) as reader:

        def build(module, writer):
            branch = writer.NewBranch("Weight", ROOT.Weight.Class())
            reader.AnalyzeWeight(branch)

        tree_reader, fin = read_tree_branch(load_delphes, tmp_path, build)
        branch = tree_reader.UseBranch("Weight")
        tree_reader.ReadEntry(0)
        assert branch.GetEntries() == 1
        assert branch.At(0).Weight == pytest.approx(HEPMCE2_EVENT_1["weight"], rel=1e-5)
        fin.Close()


def test_hepmc3_event_metadata(load_delphes, tmp_path):
    with make_reader("HepMC3", "test.hepmc3", load_delphes) as reader:

        def build(module, writer):
            branch = writer.NewBranch("Event", ROOT.HepMCEvent.Class())
            reader.AnalyzeEvent(branch, 0)

        tree_reader, fin = read_tree_branch(load_delphes, tmp_path, build)
        branch = tree_reader.UseBranch("Event")
        tree_reader.ReadEntry(0)
        assert branch.GetEntries() == 1
        check_event_fields(branch.At(0), HEPMCE3_EVENT_1)
        fin.Close()


def test_hepmc3_analyze_weight(load_delphes, tmp_path):
    with make_reader("HepMC3", "test.hepmc3", load_delphes) as reader:

        def build(module, writer):
            branch = writer.NewBranch("Weight", ROOT.Weight.Class())
            reader.AnalyzeWeight(branch)

        tree_reader, fin = read_tree_branch(load_delphes, tmp_path, build)
        branch = tree_reader.UseBranch("Weight")
        tree_reader.ReadEntry(0)
        assert branch.GetEntries() == 1
        assert branch.At(0).Weight == pytest.approx(HEPMCE3_EVENT_1["weight"], rel=1e-5)
        fin.Close()


def test_hepmc2_hepmc3_metadata_consistency(load_delphes, tmp_path):
    with (
        make_reader("HepMC2", "test.hepmc2", load_delphes) as reader2,
        make_reader("HepMC3", "test.hepmc3", load_delphes) as reader3,
    ):

        def build(module, writer):
            branch = writer.NewBranch("Event2", ROOT.HepMCEvent.Class())
            reader2.AnalyzeEvent(branch, 0)
            branch = writer.NewBranch("Event3", ROOT.HepMCEvent.Class())
            reader3.AnalyzeEvent(branch, 0)

        tree_reader, fin = read_tree_branch(load_delphes, tmp_path, build)
        e2 = tree_reader.UseBranch("Event2")
        e3 = tree_reader.UseBranch("Event3")
        tree_reader.ReadEntry(0)
        a, b = e2.At(0), e3.At(0)
        assert a.Number == b.Number
        assert a.ProcessID == b.ProcessID
        assert a.MPI == b.MPI
        assert a.Weight == pytest.approx(b.Weight, rel=1e-5)
        assert a.Scale == pytest.approx(b.Scale, rel=1e-5)
        assert a.AlphaQED == pytest.approx(b.AlphaQED, rel=1e-5)
        assert a.AlphaQCD == pytest.approx(b.AlphaQCD, rel=1e-5)
        assert a.ID1 == b.ID1
        assert a.ID2 == b.ID2
        assert a.X1 == pytest.approx(b.X1, rel=1e-5)
        assert a.X2 == pytest.approx(b.X2, rel=1e-5)
        assert a.ScalePDF == pytest.approx(b.ScalePDF, rel=1e-5)
        assert a.PDF1 == pytest.approx(b.PDF1, abs=1e-6)
        assert a.PDF2 == pytest.approx(b.PDF2, abs=1e-6)
        assert a.CrossSection == pytest.approx(b.CrossSection, rel=1e-5)
        fin.Close()


def minimal_hepmc2(momentum_unit, position_unit, event_lines=None):
    if event_lines is None:
        event_lines = [
            "E 1 1 100.0 0.116 0.0075 1 0 1 1 2 0 1 2.5",
            "U GEV MM",
            "C 10.0 1.0",
            "F 21 21 0.1 0.2 100.0 0 0",
            "V -1 0 1.0 2.0 3.0 4.0 0 2 0",
            "P 1 21 100000.0 0.0 0.0 100000.0 0 -1 0 0 -9",
            "P 2 21 0.0 0.0 -200000.0 200000.0 0 -1 0 0 -9",
        ]
        event_lines[1] = f"U {momentum_unit} {position_unit}"
    lines = [
        "HepMC::Version 2.07.00",
        "HepMC::IO_GenEvent-START_EVENT_LISTING",
    ]
    lines += event_lines
    return "\n".join(lines) + "\n"


def read_hepmc2_file(tmp_path, load_delphes, content):
    path = tmp_path / "minimal.hepmc2"
    path.write_text(content)
    with reader_context(load_delphes, "HepMC2", path) as (_, factory, _, arrays, reader):
        ok = reader.ReadEvent(factory, *arrays)
    return ok, arrays


@pytest.mark.parametrize(
    "momentum_unit,position_unit,px_scale,pos_scale",
    [
        ("GEV", "MM", 1.0, 1.0),
        ("MEV", "CM", 0.001, 10.0),
    ],
)
def test_hepmc2_units_coefficients(load_delphes, tmp_path, momentum_unit, position_unit, px_scale, pos_scale):
    content = minimal_hepmc2(momentum_unit, position_unit)
    ok, arrays = read_hepmc2_file(tmp_path, load_delphes, content)
    assert ok
    all_particles, stable, partons = arrays
    assert all_particles.GetEntries() == 2
    px = all_particles.At(0).Momentum.Px()
    pz = all_particles.At(1).Momentum.Pz()
    assert px == pytest.approx(100000.0 * px_scale, rel=1e-12)
    assert pz == pytest.approx(-200000.0 * px_scale, rel=1e-12)
    for i in range(2):
        p = all_particles.At(i).Position
        assert p.X() == pytest.approx(1.0 * pos_scale, rel=1e-12)
        assert p.Y() == pytest.approx(2.0 * pos_scale, rel=1e-12)
        assert p.Z() == pytest.approx(3.0 * pos_scale, rel=1e-12)
        assert p.T() == pytest.approx(4.0 * pos_scale, rel=1e-12)

    assert stable.GetEntries() == 0
    assert partons.GetEntries() == 2
    assert all(partons.At(i).PID == 21 and partons.At(i).Status == -1 for i in range(2))


def test_hepmc2_unknown_units_silently_unit_coefficient(load_delphes, tmp_path):
    content = minimal_hepmc2("FOO", "BAR")
    ok, arrays = read_hepmc2_file(tmp_path, load_delphes, content)
    assert ok
    assert arrays[0].GetEntries() == 2
    assert arrays[0].At(0).Momentum.Px() == pytest.approx(100000.0, rel=1e-12)
    assert arrays[0].At(0).Position.X() == pytest.approx(1.0, rel=1e-12)


def test_hepmc2_malformed_units_line_rejected(load_delphes, tmp_path):
    lines = minimal_hepmc2("GEV", "MM").splitlines()
    lines = [line for line in lines if not line.startswith("U ")]
    lines.insert(3, "U FOO")
    ok, _ = read_hepmc2_file(tmp_path, load_delphes, "\n".join(lines) + "\n")
    assert not ok


def test_hepmc2_truncated_event_line_rejected(load_delphes, tmp_path):
    lines = minimal_hepmc2("GEV", "MM").splitlines()
    lines = [line for line in lines if not line.startswith("E ")]
    lines.insert(2, "E 1 1 100.0 0.116")
    ok, _ = read_hepmc2_file(tmp_path, load_delphes, "\n".join(lines) + "\n")
    assert not ok


def test_hepmc2_malformed_particle_line_rejected(load_delphes, tmp_path):
    lines = minimal_hepmc2("GEV", "MM").splitlines()
    lines = [line for line in lines if not line.startswith("P ")]
    lines.append("P 1 21 abc 0.0 0.0 100000.0 0 -1 0 0 -9")
    ok, _ = read_hepmc2_file(tmp_path, load_delphes, "\n".join(lines) + "\n")
    assert not ok


def test_hepmc2_beam_status4_without_links(load_delphes, tmp_path):
    lines = [
        "HepMC::Version 2.07.00",
        "HepMC::IO_GenEvent-START_EVENT_LISTING",
        "E 1 1 100.0 0.116 0.0075 1 0 1 1 2 0 1 2.5",
        "U GEV MM",
        "C 10.0 1.0",
        "F 21 21 0.1 0.2 100.0 0 0",
        "V -1 0 0 0 0 0 0 1 0",
        "P 1 21 0.0 0.0 1000.0 1000.0 0 4 0 0 -7",
    ]
    ok, arrays = read_hepmc2_file(tmp_path, load_delphes, "\n".join(lines) + "\n")
    assert ok
    all_particles, stable, partons = arrays
    assert all_particles.GetEntries() == 1
    c = all_particles.At(0)
    assert c.Status == 4
    assert c.PID == 21
    assert (c.M1, c.M2, c.D1, c.D2) == (-1, -1, -1, -1)
    assert c.Momentum.Pz() == pytest.approx(1000.0, rel=1e-12)
    assert stable.GetEntries() == 0
    assert partons.GetEntries() == 1
    assert partons.At(0).PID == 21


def test_hepmc2_multi_vertex_event(load_delphes, tmp_path):
    def e(px, py, pz, m):
        return math.sqrt(px * px + py * py + pz * pz + m * m)

    lines = [
        "HepMC::Version 2.07.00",
        "HepMC::IO_GenEvent-START_EVENT_LISTING",
        "E 1 2 100.0 0.116 0.0075 1 0 5 1 2 0 1 2.5",
        "U GEV MM",
        "C 10.0 1.0",
        "F 21 21 0.1 0.2 100.0 0 0",
        "V -1 0 0.0 0.0 0.0 0.0 0 1 0",
        f"P 1 21 0.0 0.0 600.0 {e(0, 0, 600, 0):.17g} 0 -1 0 0 -3",
        "V -2 0 0.0 0.0 0.0 0.0 0 1 0",
        f"P 2 21 0.0 0.0 -400.0 {e(0, 0, 400, 0):.17g} 0 -1 0 0 -3",
        "V -3 0 0.0 0.0 0.0 0.0 0 2 0",
        f"P 3 6 10.0 0.0 500.0 {e(10, 0, 500, 173.0):.17g} 173.0 1 0 0 -5",
        f"P 4 -6 -10.0 0.0 -100.0 {e(-10, 0, -100, 173.0):.17g} 173.0 1 0 0 -6",
        "V -5 0 0.0 0.0 100.0 0.0 0 1 0",
        f"P 5 5 10.0 0.0 500.0 {e(10, 0, 500, 4.18):.17g} 4.18 1 0 0 -8",
        "V -6 0 0.0 0.0 -200.0 0.0 0 1 0",
        f"P 6 -5 -10.0 0.0 -100.0 {e(-10, 0, -100, 4.18):.17g} 4.18 1 0 0 -9",
    ]
    ok, arrays = read_hepmc2_file(tmp_path, load_delphes, "\n".join(lines) + "\n")
    assert ok
    all_particles, stable, partons = arrays
    assert all_particles.GetEntries() == 6
    p = [all_particles.At(i) for i in range(6)]
    assert [c.PID for c in p] == [21, 21, 6, -6, 5, -5]

    for i in (0, 1):
        assert (p[i].M1, p[i].M2) == (-1, -1)
        assert (p[i].D1, p[i].D2) == (2, 3)
        pos = p[i].Position
        assert (pos.X(), pos.Y(), pos.Z(), pos.T()) == (0.0, 0.0, 0.0, 0.0)

    assert (p[2].M1, p[2].M2) == (0, 1)
    assert (p[2].D1, p[2].D2) == (4, 4)
    dpos = p[2].DecayPosition
    assert (dpos.X(), dpos.Y(), dpos.Z(), dpos.T()) == (0.0, 0.0, 100.0, 0.0)
    assert (p[3].M1, p[3].M2) == (0, 1)
    assert (p[3].D1, p[3].D2) == (5, 5)
    dpos = p[3].DecayPosition
    assert (dpos.X(), dpos.Y(), dpos.Z(), dpos.T()) == (0.0, 0.0, -200.0, 0.0)

    assert (p[4].M1, p[4].M2) == (2, -1)
    assert (p[4].D1, p[4].D2) == (-1, -1)
    pos = p[4].Position
    assert (pos.X(), pos.Y(), pos.Z(), pos.T()) == (0.0, 0.0, 100.0, 0.0)
    assert (p[5].M1, p[5].M2) == (3, -1)
    assert (p[5].D1, p[5].D2) == (-1, -1)
    pos = p[5].Position
    assert (pos.X(), pos.Y(), pos.Z(), pos.T()) == (0.0, 0.0, -200.0, 0.0)

    assert stable.GetEntries() == 4
    assert partons.GetEntries() == 2
