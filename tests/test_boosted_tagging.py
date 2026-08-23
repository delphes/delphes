import math

from conftest import assert_deterministic, make_candidate, make_config, make_particle


def boosted_tagging_config(efficiency_formula, **extra):
    return make_config(
        "BoostedTagging",
        JetInputArray="Delphes/inputJets",
        ParticleInputArray="Delphes/inputParticles",
        EfficiencyFormula=efficiency_formula,
        **extra,
    )


def make_jet_candidate(factory, pt, eta, phi=0.0, soft_drop_mass=None):
    c = make_candidate(factory, pt, eta, phi, energy=pt * math.cosh(eta))
    if soft_drop_mass is not None:
        c.SoftDroppedJet.SetPtEtaPhiM(pt, eta, phi, soft_drop_mass)
    return c


def run_boosted_tagging(load_delphes, config, jets, resonances):
    module, factory, _ = load_delphes(config)
    jet_array = module.ExportArray("inputJets")
    for pt, eta, phi, mass in jets:
        jet_array.Add(make_jet_candidate(factory, pt, eta, phi, mass))
    particle_array = module.ExportArray("inputParticles")
    for pt, eta, pid in resonances:
        particle_array.Add(make_particle(factory, pt, eta, pid))
    module.Init()
    module.Process()
    return [jet_array.At(i).BoostedTag for i in range(jet_array.GetEntries())]


def test_matched_w_resonance_tagged(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0])
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 0.2, 0.0, None)], [(180.0, 0.19, 24)])
    assert tags == [1]


def test_matched_zero_efficiency_not_tagged(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 0.0])
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 0.2, 0.0, None)], [(180.0, 0.19, 24)])
    assert tags == [0]


def test_unmatched_jet_not_tagged(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0])
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 2.0, 0.0, None)], [(180.0, 1.5, 24)])
    assert tags == [0]


def test_default_efficiency_formula_fallback(load_delphes):
    tags = run_boosted_tagging(
        load_delphes, boosted_tagging_config([0, 1.0, 24, 0.0]), [(200.0, 0.2, 0.0, None)], [(180.0, 0.19, 24)]
    )
    assert tags == [0]

    tags = run_boosted_tagging(
        load_delphes, boosted_tagging_config([0, 1.0]), [(200.0, 0.2, 0.0, None)], [(180.0, 0.19, 24)]
    )
    assert tags == [1]


def test_efficiency_keyed_by_pdg(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0])
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 0.2, 0.0, None)], [(180.0, 0.19, 23)])
    assert tags == [0]


def test_soft_drop_mass_window(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0], SoftDropMassMin=50.0, SoftDropMassMax=100.0)
    jets = [(200.0, 0.2, 0.0, 80.0), (200.0, 0.2, 0.0, 20.0), (200.0, 0.2, 0.0, 120.0)]
    tags = run_boosted_tagging(load_delphes, config, jets, [(180.0, 0.19, 24)])
    assert tags == [1, 0, 0]


def test_jet_pt_min(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0], JetPTMin=100.0)
    tags = run_boosted_tagging(load_delphes, config, [(50.0, 0.2, 0.0, None)], [(180.0, 0.19, 24)])
    assert tags == [0]


def test_bit_number(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0], BitNumber=1)
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 0.2, 0.0, None)], [(180.0, 0.19, 24)])
    assert tags == [2]


def test_resonance_acceptance_cuts(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0], ResonancePTMin=100.0)
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 0.2, 0.0, None)], [(50.0, 0.19, 24)])
    assert tags == [0]

    config = boosted_tagging_config([0, 0.0, 24, 1.0], ResonanceEtaMax=1.0)
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 2.0, 0.0, None)], [(180.0, 1.9, 24)])
    assert tags == [0]


def test_empty_jet_input_no_crash(load_delphes):
    config = boosted_tagging_config([0, 1.0])
    tags = run_boosted_tagging(load_delphes, config, [], [(180.0, 0.19, 24)])
    assert tags == []


def test_soft_drop_mass_at_bounds_inclusive(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0], SoftDropMassMin=50.0, SoftDropMassMax=100.0)
    jets = [(200.0, 0.2, 0.0, 50.0), (200.0, 0.2, 0.0, 100.0)]
    tags = run_boosted_tagging(load_delphes, config, jets, [(180.0, 0.19, 24)])
    assert tags == [1, 1]


def test_jet_pt_at_min_passes(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0], JetPTMin=200.0)
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 0.2, 0.0, None)], [(180.0, 0.19, 24)])
    assert tags == [1]


def test_resonance_pt_at_min_rejected(load_delphes):
    config = boosted_tagging_config([0, 0.0, 24, 1.0], ResonancePTMin=180.0)
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 0.2, 0.0, None)], [(180.0, 0.19, 24)])
    assert tags == [0]


def test_nearest_resonance_wins(load_delphes):
    config = boosted_tagging_config([0, 0.0, 23, 0.0, 24, 1.0], DeltaR=0.5)
    tags = run_boosted_tagging(load_delphes, config, [(200.0, 0.2, 0.0, None)], [(180.0, 0.69, 23), (180.0, 0.39, 24)])
    assert tags == [1]


def test_deterministic_with_fixed_seed(load_delphes):
    config = boosted_tagging_config([0, 0.5])
    jets = [(200.0, 0.2, 0.0, None), (200.0, 1.0, 0.0, None), (200.0, 1.8, 0.0, None)]
    resonances = [(180.0, 0.19, 24), (180.0, 0.99, 24), (180.0, 1.79, 24)]
    assert_deterministic(
        lambda: run_boosted_tagging(load_delphes, config, jets, resonances),
        extract=lambda tags: tuple(tags),
    )
