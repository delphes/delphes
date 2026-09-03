import math

import pytest
from .conftest import assert_deterministic, build_config, candidate_snapshots, make_candidate


def make_module_config(**extra):
    return build_config(
        "DenseTrackFilter",
        {
            "TrackInputArray": "Delphes/inputTracks",
            "TrackOutputArray": "tracks",
            "ChargedHadronOutputArray": "chargedHadrons",
            "ElectronOutputArray": "electrons",
            "MuonOutputArray": "muons",
            "EtaPhiRes": 0.0,
            "EtaPhiBins": [[-1.0, 0.0, 1.0, 2.0], [4]],
        },
        **extra,
    )


def run_module_test(run_generic, config, particles):
    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        for pt, eta, phi, charge, pid in particles:
            energy = pt
            r = 1000.0
            z = r / math.tan(2.0 * math.atan(math.exp(-eta)))
            x = r * math.cos(phi)
            y = r * math.sin(phi)

            c = make_candidate(factory, pt, eta, phi, energy=energy, pid=pid, charge=charge)
            c.Position.SetXYZT(x, y, z, 0.0)

            mother = make_candidate(factory, pt, eta, phi, energy=energy, pid=pid, charge=charge)
            mother.Position.SetXYZT(x, y, z, 0.0)
            c.AddCandidate(mother)

            input_array.Add(c)

    return run_generic(
        config,
        setup=setup,
        outputs={
            "tracks": "TestModule/tracks",
            "chargedHadrons": "TestModule/chargedHadrons",
            "electrons": "TestModule/electrons",
            "muons": "TestModule/muons",
        },
    )


def test_keeps_single_track(run_generic):
    results = run_module_test(run_generic, make_module_config(), [(50.0, 0.5, 0.5, 1, 211)])
    assert results["tracks"].GetEntries() == 1
    assert results["chargedHadrons"].GetEntries() == 1


def test_keeps_leading_track(run_generic):
    results = run_module_test(run_generic, make_module_config(), [(50.0, 0.5, 0.5, 1, 211), (10.0, 0.5, 0.5, 1, 211)])
    assert results["tracks"].GetEntries() == 1
    assert results["tracks"].At(0).Momentum.Pt() == pytest.approx(50.0, rel=1e-3)


def test_electron_classification(run_generic):
    results = run_module_test(run_generic, make_module_config(), [(50.0, 0.5, 0.5, -1, 11)])
    assert results["electrons"].GetEntries() == 1
    assert results["chargedHadrons"].GetEntries() == 0


def test_muon_classification(run_generic):
    results = run_module_test(run_generic, make_module_config(), [(50.0, 0.5, 0.5, -1, 13)])
    assert results["muons"].GetEntries() == 1
    assert results["chargedHadrons"].GetEntries() == 0


def test_empty_input(run_generic):
    results = run_module_test(run_generic, make_module_config(), [])
    assert results["tracks"].GetEntries() == 0
    assert results["chargedHadrons"].GetEntries() == 0


def make_track_at_position(factory, x, y, z, pt=50.0, charge=1, pid=211):
    c = make_candidate(factory, pt, 0.5, 0.0, pid=pid, charge=charge)
    c.Position.SetXYZT(x, y, z, 0.0)
    mother = make_candidate(factory, pt, 0.5, 0.0, pid=pid, charge=charge)
    mother.Position.SetXYZT(x, y, z, 0.0)
    c.AddCandidate(mother)
    return c


def run_position_test(run_generic, positions, eta_bins=None):
    if eta_bins is not None:
        config = build_config(
            "DenseTrackFilter",
            {
                "TrackInputArray": "Delphes/inputTracks",
                "TrackOutputArray": "tracks",
                "ChargedHadronOutputArray": "chargedHadrons",
                "ElectronOutputArray": "electrons",
                "MuonOutputArray": "muons",
                "EtaPhiRes": 0.0,
                "EtaPhiBins": [eta_bins, [4]],
            },
        )
    else:
        config = make_module_config()

    def setup(module, factory):
        input_array = module.ExportArray("inputTracks")
        for x, y, z in positions:
            input_array.Add(make_track_at_position(factory, x, y, z))

    return run_generic(config, setup=setup, outputs=("TestModule/tracks",))


def test_phi_on_internal_edge_kept(run_generic):
    results = run_position_test(run_generic, [(1000.0, 0.0, 0.0)])
    assert results.GetEntries() == 1


def test_phi_on_last_edge_kept_in_degenerate_bin(run_generic):
    results = run_position_test(run_generic, [(-1000.0, 0.0, 0.0)])
    assert results.GetEntries() == 1


def test_eta_on_first_edge_dropped(run_generic):
    results = run_position_test(run_generic, [(1000.0, 0.0, 0.0)], eta_bins=[0.0, 1.0, 2.0])
    assert results.GetEntries() == 0


def test_eta_on_internal_edge_kept(run_generic):
    results = run_position_test(run_generic, [(1000.0, 0.0, 0.0)])
    assert results.GetEntries() == 1


def test_eta_beyond_last_edge_dropped(run_generic):
    z = 1000.0 * math.sinh(2.5)
    results = run_position_test(run_generic, [(1000.0, 0.0, z)])
    assert results.GetEntries() == 0


def test_repeated_process_resets_state(load_delphes):
    config = make_module_config()
    module, factory, _ = load_delphes(config)
    input_array = module.ExportArray("inputTracks")
    for pt, eta, phi, charge, pid in [(50.0, 0.5, 0.5, 1, 211)]:
        c = make_candidate(factory, pt, eta, phi, pid=pid, charge=charge)
        c.AddCandidate(make_candidate(factory, pt, eta, phi, pid=pid, charge=charge))
        input_array.Add(c)
    module.Init()
    module.Process()
    tracks = module.ImportArray("TestModule/tracks")
    first = tracks.GetEntries()
    assert first == 1
    module.Process()
    assert tracks.GetEntries() == 2 * first


def test_deterministic_with_fixed_seed(run_generic):
    config = make_module_config(EtaPhiRes=0.01)
    particles = [(50.0, 0.5, 0.1, 1, 211), (40.0, 1.5, 2.0, 1, 211), (30.0, 0.5, 3.0, 1, 211)]
    assert_deterministic(
        lambda: run_module_test(run_generic, config, particles),
        extract=lambda results: candidate_snapshots(results["tracks"], ("Position",)),
    )


def test_extreme_eta_position_outside_bins_dropped(run_generic):
    for eta in (5.0, -5.0):
        results = run_module_test(run_generic, make_module_config(), [(50.0, eta, 0.5, 1, 211)])
        assert results["tracks"].GetEntries() == 0
        assert results["chargedHadrons"].GetEntries() == 0


def test_near_zero_pt_track_survives_bin(run_generic):
    results = run_module_test(run_generic, make_module_config(), [(1e-6, 0.5, 0.5, 1, 211)])
    assert results["tracks"].GetEntries() == 1
    assert results["chargedHadrons"].GetEntries() == 1
    assert results["tracks"].At(0).Momentum.Pt() == pytest.approx(1e-6, rel=1e-6)


def test_extreme_pt_wins_bin(run_generic):
    results = run_module_test(run_generic, make_module_config(), [(1e4, 0.5, 0.5, 1, 211), (50.0, 0.5, 0.5, 1, 211)])
    assert results["tracks"].GetEntries() == 1
    assert results["tracks"].At(0).Momentum.Pt() == pytest.approx(1e4, rel=1e-6)
