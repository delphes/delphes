import math

import pytest
from .conftest import assert_deterministic, build_config, candidate_snapshots, make_candidate, run_repeated

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

CALORIMETER_DEFAULTS = dict(
    BASE_DEFAULTS,
    ECalResolutionFormula=0.0,
    HCalResolutionFormula=0.0,
    ECalEnergyMin=0.0,
    HCalEnergyMin=0.0,
    ECalEnergySignificanceMin=0.0,
    HCalEnergySignificanceMin=0.0,
    TimingEnergyMin=0.0,
    EFlowNeutralHadronOutputArray="eflowNeutralHadrons",
)


CALORIMETERS = [
    ("Calorimeter", CALORIMETER_DEFAULTS),
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

        assert results.GetEntries() == 1

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
            extra = {"IsEcal": False, "EnergyFraction": [0, 1.0, 11, 0.0, 22, 0.0]}
        else:
            extra = {"EnergyFraction": [0, [0.0, 1.0], 130, [0.0, 1.0]]}
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

    def test_particle_on_phi_bin_edge(self, run_calorimeter, class_name, defaults):
        config = build_config(class_name, defaults)
        results = run_calorimeter(config, particles=[(50.0, 0.5, -math.pi / 2.0, 22)])
        assert results.GetEntries() == 1

    def test_smear_tower_center_branch(self, run_calorimeter, class_name, defaults):
        config = build_config(class_name, defaults, SmearTowerCenter=True)
        results = run_calorimeter(config, particles=[(50.0, 0.5, 0.0, 22)])
        assert results.GetEntries() == 1
        tower = results.At(0)
        assert tower.Edges[0] == pytest.approx(0.0, abs=1e-6)
        assert tower.Edges[1] == pytest.approx(1.0, abs=1e-6)


def run_calo_outputs(run_generic, config, particles, tracks=()):
    def setup(module, factory):
        particle_array = module.ExportArray("inputParticles")
        for pt, eta, phi, pid, t in particles:
            c = make_candidate(factory, pt, eta, phi, pid=pid, charge=0)
            c.Position.SetPtEtaPhiE(1.0, eta, phi, t)
            particle_array.Add(c)
        track_array = module.ExportArray("inputTracks")
        for pt, eta, phi, charge, pid, res, t in tracks:
            c = make_candidate(factory, pt, eta, phi, pid=pid, charge=charge)
            c.Position.SetPtEtaPhiE(1.0, eta, phi, t)
            c.TrackResolution = res
            track_array.Add(c)

    return run_generic(
        config,
        setup=setup,
        outputs={
            "towers": "TestModule/towers",
            "eflowTracks": "TestModule/eflowTracks",
            "eflowPhotons": "TestModule/eflowPhotons",
        },
    )


def test_calorimeter_eflow_photon_and_track(run_generic):
    config = build_config("Calorimeter", dict(CALORIMETERS[0][1], EFlowPhotonOutputArray="eflowPhotons"))
    results = run_calo_outputs(
        run_generic, config, particles=[(50.0, 0.5, 0.0, 22, 1000.0)], tracks=[(20.0, 0.5, 0.0, 1, 211, 0.01, 500.0)]
    )
    assert results["towers"].GetEntries() == 1
    assert results["eflowTracks"].GetEntries() == 1
    assert results["eflowPhotons"].GetEntries() == 1


def test_calorimeter_timing_gate(run_generic):
    low = run_calo_outputs(
        run_generic,
        build_config(
            "Calorimeter", dict(CALORIMETERS[0][1], EFlowPhotonOutputArray="eflowPhotons", TimingEnergyMin=0.0)
        ),
        particles=[(50.0, 0.5, 0.0, 22, 1000.0)],
    )
    assert low["towers"].At(0).Position.T() == pytest.approx(1000.0, rel=1e-3)
    high = run_calo_outputs(
        run_generic,
        build_config(
            "Calorimeter", dict(CALORIMETERS[0][1], EFlowPhotonOutputArray="eflowPhotons", TimingEnergyMin=1000.0)
        ),
        particles=[(50.0, 0.5, 0.0, 22, 1000.0)],
    )
    assert high["towers"].At(0).Position.T() == pytest.approx(999999.9, abs=0.1)


def test_dual_readout_log_normal_smearing(run_generic):
    config = build_config(
        "DualReadoutCalorimeter",
        dict(CALORIMETERS[3][1], SmearLogNormal=True, ECalResolutionFormula=0.1, HCalResolutionFormula=0.1),
    )
    results = run_calo_outputs(run_generic, config, particles=[(50.0, 0.5, 0.0, 22, 0.0)])
    assert results["towers"].GetEntries() == 1
    assert results["towers"].At(0).Momentum.E() != pytest.approx(50.0 * math.cosh(0.5), rel=1e-6)


@pytest.mark.parametrize(
    "class_name,defaults,extra",
    [
        ("Calorimeter", CALORIMETERS[0][1], {"ECalEnergyMin": 100.0, "HCalEnergyMin": 100.0}),
        ("SimpleCalorimeter", CALORIMETERS[2][1], {"EnergyMin": 100.0}),
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
        ("SimpleCalorimeter", CALORIMETERS[2][1], {"IsEcal": False, "EnergyFraction": [0, 1.0, 11, 0.0, 22, 0.0]}),
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
    assert results.GetEntries() == 1
    tower = results.At(0)
    assert tower.Etrk > 0


CALO_RESOLUTION_OVERRIDE = {
    "Calorimeter": {"ECalResolutionFormula": 0.1},
    "OldCalorimeter": {"ECalResolutionFormula": 0.1},
    "SimpleCalorimeter": {"ResolutionFormula": 0.1},
    "DualReadoutCalorimeter": {"ECalResolutionFormula": 0.1, "HCalResolutionFormula": 0.1},
}


def test_charged_hadron_ecal_hcal_split(run_calorimeter):
    config = build_config("Calorimeter", dict(CALORIMETERS[0][1], EnergyFraction=[211, [0.1, 0.9]]))
    pt, eta = 50.0, 0.5
    energy = pt * math.cosh(eta)
    results = run_calorimeter(config, particles=[(pt, eta, 0.0, 211)])
    assert results.GetEntries() == 1
    tower = results.At(0)

    assert tower.Eem == pytest.approx(0.1 * energy, abs=1e-5)
    assert tower.Ehad == pytest.approx(0.9 * energy, abs=1e-5)

    assert tower.Momentum.E() == pytest.approx(energy, rel=1e-9)


@pytest.mark.parametrize("class_name,defaults", CALORIMETERS, ids=CALORIMETER_IDS)
def test_multiple_particles_merge_into_single_tower(run_calorimeter, class_name, defaults):
    config = build_config(class_name, defaults)
    results = run_calorimeter(config, particles=[(30.0, 0.5, 0.0, 22), (20.0, 0.5, 0.0, 22)])
    assert results.GetEntries() == 1
    tower = results.At(0)
    e_total = 30.0 * math.cosh(0.5) + 20.0 * math.cosh(0.5)
    assert tower.Momentum.E() == pytest.approx(e_total, rel=1e-9)
    assert tower.GetCandidates().GetEntries() == 2


@pytest.mark.parametrize("class_name,defaults", CALORIMETERS, ids=CALORIMETER_IDS)
def test_eta_acceptance_edges(run_calorimeter, class_name, defaults):
    config = build_config(class_name, defaults)
    inside = run_calorimeter(config, particles=[(50.0, 1.5, 0.0, 22)])
    assert inside.GetEntries() == 1
    on_edge = run_calorimeter(config, particles=[(50.0, 1.0, 0.0, 22)])
    assert on_edge.GetEntries() == 1
    assert on_edge.At(0).Edges[0] == pytest.approx(0.0, abs=1e-6)
    assert on_edge.At(0).Edges[1] == pytest.approx(1.0, abs=1e-6)
    beyond = run_calorimeter(config, particles=[(50.0, 2.5, 0.0, 22)])
    assert beyond.GetEntries() == 0
    on_last_edge = run_calorimeter(config, particles=[(50.0, 2.0, 0.0, 22)])
    assert on_last_edge.GetEntries() == 0


def test_ecal_energy_significance_min_drops_insignificant_tower(run_generic):
    defaults = dict(CALORIMETERS[0][1], ECalResolutionFormula=1.0, ECalEnergySignificanceMin=5.0)

    def setup(module, factory):
        particle_array = module.ExportArray("inputParticles")
        for pt, eta in ((2.0, 0.5), (10.0, 1.5)):
            c = make_candidate(factory, pt, eta, 0.0, pid=22, charge=0)
            c.Position.SetPtEtaPhiE(1.0, eta, 0.0, 0)
            particle_array.Add(c)
        module.ExportArray("inputTracks")

    config = build_config("Calorimeter", defaults)
    results = run_generic(config, setup=setup, outputs=("TestModule/towers",))
    assert results.GetEntries() == 1
    tower = results.At(0)

    assert tower.Edges[0] == pytest.approx(1.0, abs=1e-6)
    assert 20.0 < tower.Momentum.E() < 27.0

    control = build_config("Calorimeter", dict(CALORIMETERS[0][1]))
    results = run_generic(control, setup=setup, outputs=("TestModule/towers",))
    assert results.GetEntries() == 2


def dual_readout_gaussian_config(formula):
    return build_config(
        "DualReadoutCalorimeter",
        dict(CALORIMETERS[3][1], SmearLogNormal=False, ECalResolutionFormula=formula),
    )


def truncated_normal_moments(mu, sigma):
    a = mu / sigma
    phi = math.exp(-0.5 * a * a) / math.sqrt(2.0 * math.pi)
    cdf = 0.5 * (1.0 + math.erf(a / math.sqrt(2.0)))
    mean = mu + sigma * phi / cdf
    var = sigma * sigma * (1.0 + a * phi / cdf - phi * phi / (cdf * cdf))
    return mean, var


def test_dual_readout_gaussian_path(run_generic, load_delphes):
    pt, eta = 50.0, 0.0

    def setup(module, factory):
        particle_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, pt, eta, 0.0, pid=22, charge=0)
        c.Position.SetPtEtaPhiE(1.0, eta, 0.0, 0)
        particle_array.Add(c)
        module.ExportArray("inputTracks")

    config = dual_readout_gaussian_config("0.5*energy")

    e1 = run_generic(config, setup=setup, outputs=("TestModule/towers",)).At(0).Momentum.E()
    e2 = run_generic(config, setup=setup, outputs=("TestModule/towers",)).At(0).Momentum.E()
    assert e1 == e2

    n = 500
    towers_list = run_repeated(load_delphes, config, n, setup, output="TestModule/towers")
    energies = [towers.At(0).Momentum.E() for towers in towers_list if towers.GetEntries() == 1]

    assert len(energies) < n
    assert len(energies) > 0.9 * n
    exp_mean, exp_var = truncated_normal_moments(pt, 0.5 * pt)
    m = sum(energies) / len(energies)
    var = sum((e - m) ** 2 for e in energies) / (len(energies) - 1)
    assert m == pytest.approx(exp_mean, abs=4.0)
    assert math.sqrt(var) == pytest.approx(math.sqrt(exp_var), abs=2.5)

    log_config = build_config(
        "DualReadoutCalorimeter",
        dict(CALORIMETERS[3][1], SmearLogNormal=True, ECalResolutionFormula="0.5*energy"),
    )
    towers_list = run_repeated(load_delphes, log_config, n, setup, output="TestModule/towers")
    kept = sum(towers.GetEntries() for towers in towers_list)
    assert kept == n


@pytest.mark.parametrize("class_name,defaults", CALORIMETERS, ids=CALORIMETER_IDS)
def test_deterministic_with_fixed_seed(run_calorimeter, class_name, defaults):
    config = build_config(class_name, dict(defaults, **CALO_RESOLUTION_OVERRIDE[class_name]))
    particles = [(50.0, 0.5, 0.0, 22), (60.0, 1.5, 0.2, 22), (40.0, -0.5, 0.4, 22)]
    assert_deterministic(
        lambda: run_calorimeter(config, particles),
        extract=lambda towers: candidate_snapshots(towers, ("Momentum", "Eem", "Ehad")),
    )
