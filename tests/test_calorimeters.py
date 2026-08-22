import math
import pytest
from conftest import build_config

BASE_DEFAULTS = {
    "ParticleInputArray": "Delphes/inputParticles",
    "TrackInputArray": "Delphes/inputTracks",
    "TowerOutputArray": "towers",
    "PhotonOutputArray": "photons",
    "EFlowTrackOutputArray": "eflowTracks",
    "EFlowTowerOutputArray": "eflowTowers",
    "SmearTowerCenter": False,
    "EtaPhiBins": [[-1.0, 0.0, 1.0, 2.0], [4]],
    "EnergyFraction": [0, [0.0, 1.0], 11, [1.0, 0.0], 22, [1.0, 0.0]],
}

CALORIMETERS = [
    (
        "Calorimeter",
        dict(
            BASE_DEFAULTS,
            ECalResolutionFormula=0.0,
            HCalResolutionFormula=0.0,
            ECalEnergyMin=0.0,
            HCalEnergyMin=0.0,
            ECalEnergySignificanceMin=0.0,
            HCalEnergySignificanceMin=0.0,
            TimingEnergyMin=0.0,
            EFlowNeutralHadronOutputArray="eflowNeutralHadrons",
        ),
    ),
    (
        "OldCalorimeter",
        dict(
            BASE_DEFAULTS,
            ECalResolutionFormula=0.0,
            HCalResolutionFormula=0.0,
        ),
    ),
    (
        "SimpleCalorimeter",
        dict(
            BASE_DEFAULTS,
            ResolutionFormula=0.0,
            EnergyMin=0.0,
            EnergySignificanceMin=0.0,
            IsEcal=True,
            EnergyFraction=[0, 0.0, 11, 1.0, 22, 1.0],
        ),
    ),
    (
        "DualReadoutCalorimeter",
        dict(
            BASE_DEFAULTS,
            SmearLogNormal=False,
            ECalResolutionFormula=0.0,
            HCalResolutionFormula=0.0,
            ECalMinSignificance=0.0,
            HCalMinSignificance=0.0,
            TimingEnergyMin=0.0,
            EFlowPhotonOutputArray="eflowPhotons",
            EFlowNeutralHadronOutputArray="eflowNeutralHadrons",
        ),
    ),
]

CALORIMETER_IDS = [name for name, _ in CALORIMETERS]


@pytest.mark.parametrize("class_name,defaults", CALORIMETERS, ids=CALORIMETER_IDS)
class TestCalorimeters:
    def test_photon_creates_tower(self, run_calorimeter, class_name, defaults):
        config = build_config(class_name, defaults)
        results = run_calorimeter(config, particles=[(50.0, 0.5, 0.0, 22)])
        assert results.GetEntries() >= 1

    def test_energy_conservation(self, run_calorimeter, class_name, defaults):
        config = build_config(class_name, defaults)

        pt = 50.0
        eta = 0.5
        expected_energy = pt * math.cosh(eta)
        results = run_calorimeter(config, particles=[(pt, eta, 0.0, 22)])
        total_energy = 0.0
        for i in range(results.GetEntries()):
            total_energy += results.At(i).Momentum.E()
        assert total_energy == pytest.approx(expected_energy, rel=1e-3)

    def test_ecal_sets_eem(self, run_calorimeter, class_name, defaults):
        config = build_config(class_name, defaults)
        results = run_calorimeter(config, particles=[(50.0, 0.5, 0.0, 22)])
        tower = results.At(0)
        assert tower.Eem > 0
        assert tower.Ehad == 0

    def test_hcal_sets_ehad(self, run_calorimeter, class_name, defaults):
        if class_name == "SimpleCalorimeter":
            extra = dict(IsEcal=False, EnergyFraction=[0, 1.0, 11, 0.0, 22, 0.0])
        else:
            extra = dict(EnergyFraction=[0, [0.0, 1.0], 130, [0.0, 1.0]])
        config = build_config(class_name, defaults, **extra)
        results = run_calorimeter(config, particles=[(50.0, 0.5, 0.0, 130)])
        tower = results.At(0)
        assert tower.Ehad > 0
        assert tower.Eem == 0

    def test_tower_eta_phi_edges(self, run_calorimeter, class_name, defaults):
        config = build_config(class_name, defaults)
        results = run_calorimeter(config, particles=[(50.0, 1.0, 0.0, 22)])
        tower = results.At(0)
        assert tower.Edges[0] == pytest.approx(0.0, abs=1e-6)
        assert tower.Edges[1] == pytest.approx(1.0, abs=1e-6)
        assert tower.Edges[2] == pytest.approx(-1.570796, abs=1e-6)
        assert tower.Edges[3] == pytest.approx(0.0, abs=1e-6)


@pytest.mark.parametrize(
    "class_name,defaults,extra",
    [
        ("Calorimeter", CALORIMETERS[0][1], dict(ECalEnergyMin=100.0, HCalEnergyMin=100.0)),
        ("SimpleCalorimeter", CALORIMETERS[2][1], dict(EnergyMin=100.0)),
    ],
    ids=["Calorimeter", "SimpleCalorimeter"],
)
def test_energy_min_filter(run_calorimeter, class_name, defaults, extra):
    config = build_config(class_name, defaults, **extra)
    results = run_calorimeter(config, particles=[(50.0, 0.5, 0.0, 22)])
    assert results.GetEntries() == 0


@pytest.mark.parametrize(
    "class_name,defaults,extra",
    [
        ("SimpleCalorimeter", CALORIMETERS[2][1], dict(IsEcal=False, EnergyFraction=[0, 1.0, 11, 0.0, 22, 0.0])),
        ("DualReadoutCalorimeter", CALORIMETERS[3][1], {}),
    ],
    ids=["SimpleCalorimeter", "DualReadoutCalorimeter"],
)
def test_charged_track_in_tower(run_calorimeter, class_name, defaults, extra):
    config = build_config(class_name, defaults, **extra)
    results = run_calorimeter(
        config,
        particles=[(50.0, 0.5, 0.0, 130)],
        tracks=[(50.0, 0.5, 0.0, 1, 211, 0.01)],
    )
    assert results.GetEntries() >= 1
    tower = results.At(0)
    assert tower.Etrk > 0
