import pytest
from .conftest import (
    DATA_DIR,
    assert_deterministic,
    candidate_snapshots,
    make_config,
    pythia8_pileup_available,
    run_pileup,
)

CMND_FILE = str(DATA_DIR / "test.cmnd")

PION = [(50.0, 1)]


def pileup_config(**extra):
    params = {
        "ConfigFile": CMND_FILE,
        "PileUpDistribution": 2,
        "VertexD": "1.0",
    }
    params.update(extra)
    return make_config("PileUpMergerPythia8", **params)


@pytest.mark.skipif(not pythia8_pileup_available, reason="PileUpMergerPythia8 not built (no Pythia8)")
def test_output_arrays(load_delphes):
    parts, verts = run_pileup(load_delphes, pileup_config(MeanPileUp=1.0), PION)

    assert verts.GetEntries() == 2
    assert verts.At(0).IsPU == 0
    assert verts.At(1).IsPU == 1

    assert parts.GetEntries() > 1
    assert parts.At(0).IsPU == 0
    assert parts.At(0).PID == 211

    assert any(parts.At(i).IsPU == 1 for i in range(1, parts.GetEntries()))


@pytest.mark.skipif(not pythia8_pileup_available, reason="PileUpMergerPythia8 not built (no Pythia8)")
def test_pt_min_cut_removes_soft_pileup(load_delphes):
    parts_lo, _ = run_pileup(load_delphes, pileup_config(MeanPileUp=1.0, PTMin=0.0), PION)
    parts_hi, _ = run_pileup(load_delphes, pileup_config(MeanPileUp=1.0, PTMin=50.0), PION)
    n_lo = sum(1 for i in range(1, parts_lo.GetEntries()) if parts_lo.At(i).IsPU == 1)
    n_hi = sum(1 for i in range(1, parts_hi.GetEntries()) if parts_hi.At(i).IsPU == 1)
    assert n_lo > 0
    assert n_hi < n_lo


@pytest.mark.skipif(not pythia8_pileup_available, reason="PileUpMergerPythia8 not built (no Pythia8)")
def test_deterministic_with_fixed_seed(load_delphes):
    config = pileup_config(MeanPileUp=1.0)

    assert_deterministic(
        lambda: run_pileup(load_delphes, config, PION),
        extract=lambda pair: tuple(candidate_snapshots(arr, ("PID", "IsPU", "Momentum", "Position")) for arr in pair),
        abs_tol=1e-6,
    )
