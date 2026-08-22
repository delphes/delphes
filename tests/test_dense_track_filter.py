import math
import pytest
from conftest import build_config, make_candidate


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
