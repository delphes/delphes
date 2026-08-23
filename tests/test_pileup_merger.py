import math

import pytest
from conftest import (
    C_LIGHT,
    MINBIAS_FILE,
    assert_deterministic,
    candidate_snapshots,
    make_candidate,
    make_config,
    run_pileup,
)

Z_SPREAD_MM = 0.15 * 1.0e3
T_SPREAD_MM = 1.5e-9 * C_LIGHT * 1.0e3


def pileup_config(**extra):
    params = {
        "PileUpFile": MINBIAS_FILE,
        "PileUpDistribution": 2,
        "VertexDistributionFormula": "1.0",
    }
    params.update(extra)
    return make_config("PileUpMerger", **params)


SNAPSHOT_FIELDS = ("PID", "IsPU", "Momentum", "Position")


INPUT_PARTICLES = [(50.0, 1), (30.0, -1)]


def test_zero_mean_pileup_no_pileup_particles(load_delphes):
    parts, verts = run_pileup(load_delphes, pileup_config(MeanPileUp=0.0), INPUT_PARTICLES)
    assert parts.GetEntries() == 2
    assert all(parts.At(i).IsPU == 0 for i in range(parts.GetEntries()))
    assert verts.GetEntries() == 1
    vertex = verts.At(0)
    assert vertex.IsPU == 0
    assert vertex.ClusterNDF == 2
    assert vertex.SumPT2 == pytest.approx(50.0**2 + 30.0**2, rel=1e-6)


def test_positive_mean_pileup_adds_pileup_particles(load_delphes):
    parts, verts = run_pileup(load_delphes, pileup_config(MeanPileUp=2.0), INPUT_PARTICLES)
    assert verts.GetEntries() == 3
    assert parts.GetEntries() > 2

    assert all(parts.At(i).IsPU == 0 for i in range(2))

    assert any(parts.At(i).IsPU == 1 for i in range(2, parts.GetEntries()))
    assert all(verts.At(i).IsPU == 1 for i in range(1, 3))


def test_primary_vertex_position_from_formula(load_delphes):
    parts, verts = run_pileup(load_delphes, pileup_config(MeanPileUp=0.0), INPUT_PARTICLES)
    vertex = verts.At(0)

    assert abs(vertex.Position.Z()) <= Z_SPREAD_MM
    assert abs(vertex.Position.T()) <= T_SPREAD_MM

    particle = parts.At(0)
    assert particle.Position.Z() == pytest.approx(vertex.Position.Z(), rel=1e-6)
    assert particle.Position.T() == pytest.approx(vertex.Position.T(), rel=1e-6)


def test_output_beam_spot_offset(load_delphes):
    parts_a, _ = run_pileup(load_delphes, pileup_config(MeanPileUp=2.0, OutputBeamSpotX=0.0), INPUT_PARTICLES)
    parts_b, _ = run_pileup(load_delphes, pileup_config(MeanPileUp=2.0, OutputBeamSpotX=2.0), INPUT_PARTICLES)
    snap_a = candidate_snapshots(parts_a, SNAPSHOT_FIELDS)
    snap_b = candidate_snapshots(parts_b, SNAPSHOT_FIELDS)
    assert len(snap_a) == len(snap_b)

    for i in range(2):
        assert snap_b[i][6] == pytest.approx(snap_a[i][6], abs=1e-9)

    for i in range(2, len(snap_a)):
        assert snap_b[i][6] - snap_a[i][6] == pytest.approx(2.0, abs=1e-6)
        assert snap_b[i][7] == pytest.approx(snap_a[i][7], abs=1e-6)


def test_input_beam_spot_offset(load_delphes):
    parts_a, _ = run_pileup(load_delphes, pileup_config(MeanPileUp=2.0, InputBeamSpotX=0.0), INPUT_PARTICLES)
    parts_b, _ = run_pileup(load_delphes, pileup_config(MeanPileUp=2.0, InputBeamSpotX=2.0), INPUT_PARTICLES)
    snap_a = candidate_snapshots(parts_a, SNAPSHOT_FIELDS)
    snap_b = candidate_snapshots(parts_b, SNAPSHOT_FIELDS)

    for i in range(2, len(snap_a)):
        dx = snap_b[i][6] - snap_a[i][6]
        dy = snap_b[i][7] - snap_a[i][7]
        assert math.hypot(dx, dy) == pytest.approx(2.0, abs=1e-6)
    for i in range(2):
        assert snap_b[i][6] == pytest.approx(snap_a[i][6], abs=1e-9)


def test_deterministic_with_fixed_seed(load_delphes):
    config = pileup_config(MeanPileUp=2.0)
    assert_deterministic(
        lambda: run_pileup(load_delphes, config, INPUT_PARTICLES),
        extract=lambda pair: tuple(candidate_snapshots(arr, SNAPSHOT_FIELDS) for arr in pair),
    )


def test_empty_input(load_delphes):
    parts, verts = run_pileup(load_delphes, pileup_config(MeanPileUp=0.0), [])
    assert parts.GetEntries() == 0
    assert verts.GetEntries() == 1
    assert verts.At(0).ClusterNDF == 0


def test_repeated_process_resets_state(load_delphes):
    config = pileup_config(MeanPileUp=1.0)
    module, factory, _ = load_delphes(config)
    input_array = module.ExportArray("inputParticles")
    for pt, charge in INPUT_PARTICLES:
        c = make_candidate(factory, pt, 0.5, pid=211, charge=charge)
        c.Position.SetXYZT(0.0, 0.0, 0.0, 0.0)
        input_array.Add(c)
    module.Init()
    module.Process()
    parts = module.ImportArray("TestModule/stableParticles")
    verts = module.ImportArray("TestModule/vertices")
    first_parts = parts.GetEntries()
    first_verts = verts.GetEntries()
    assert first_verts == 2
    module.Process()
    assert verts.GetEntries() == 2 * first_verts
    assert parts.GetEntries() > first_parts
