import pytest
from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "LeptonDressing",
        {
            "CandidateInputArray": "Delphes/inputCandidates",
            "DressingInputArray": "Delphes/inputDressing",
            "OutputArray": "outputParticles",
            "DeltaRMax": 0.5,
        },
        **extra,
    )


def run_dressing_test(run_generic, config, candidate_pt=50.0, dressing_pt=5.0, delta_r=0.1):
    def setup(module, factory):
        candidate_array = module.ExportArray("inputCandidates")
        c = make_candidate(factory, candidate_pt, 0.5)
        candidate_array.Add(c)

        dressing_array = module.ExportArray("inputDressing")
        d = make_candidate(factory, dressing_pt, 0.5, 0.0 + delta_r)
        dressing_array.Add(d)

    return run_generic(config, setup=setup, outputs=("TestModule/outputParticles",))


def test_dressing_adds_momentum(run_generic):
    output = run_dressing_test(run_generic, make_module_config(), candidate_pt=50.0, dressing_pt=5.0, delta_r=0.1)
    assert output.GetEntries() == 1
    dressed = output.At(0)
    assert dressed.Momentum.Pt() > 50.0


def test_dressing_no_match(run_generic):
    output = run_dressing_test(
        run_generic, make_module_config(DeltaRMax=0.1), candidate_pt=50.0, dressing_pt=5.0, delta_r=3.0
    )
    dressed = output.At(0)
    assert dressed.Momentum.Pt() == pytest.approx(50.0, rel=1e-3)


def test_dressing_preserves_candidate(run_generic):
    output = run_dressing_test(run_generic, make_module_config(), candidate_pt=50.0, dressing_pt=5.0, delta_r=0.1)
    dressed = output.At(0)
    assert dressed.Momentum.Eta() == pytest.approx(0.5, rel=1e-3)


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputCandidates")
        module.ExportArray("inputDressing")

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputParticles",))
    assert output.GetEntries() == 0


def run_with_refs(load_delphes, config, dressing_specs):
    module, factory, _ = load_delphes(config)
    candidate_array = module.ExportArray("inputCandidates")
    c = make_candidate(factory, 50.0, 0.5)
    candidate_array.Add(c)
    dressing_array = module.ExportArray("inputDressing")
    dressings = []
    for pt, eta, phi in dressing_specs:
        d = make_candidate(factory, pt, eta, phi)
        dressing_array.Add(d)
        dressings.append(d)
    module.Init()
    module.Process()
    out = module.ImportArray("TestModule/outputParticles")
    return c, dressings, out


def test_dressed_four_momentum_is_exact_sum(load_delphes):
    c, dressings, out = run_with_refs(
        load_delphes, make_module_config(DeltaRMax=0.5), [(5.0, 0.5, 0.1), (3.0, 0.6, -0.2)]
    )
    assert out.GetEntries() == 1
    dressed = out.At(0).Momentum
    expected = c.Momentum
    for d in dressings:
        expected += d.Momentum
    assert dressed.Px() == pytest.approx(expected.Px(), rel=1e-12)
    assert dressed.Py() == pytest.approx(expected.Py(), rel=1e-12)
    assert dressed.Pz() == pytest.approx(expected.Pz(), rel=1e-12)
    assert dressed.E() == pytest.approx(expected.E(), rel=1e-12)


def test_lepton_without_soft_particles_unchanged(load_delphes):
    c, _, out = run_with_refs(load_delphes, make_module_config(DeltaRMax=0.5), [(5.0, 3.0, 0.1)])
    assert out.At(0).Momentum.E() == c.Momentum.E()
    assert out.At(0).Momentum.Pt() == c.Momentum.Pt()


def test_dressing_just_inside_delta_r_included(load_delphes):
    c, dressings, out = run_with_refs(load_delphes, make_module_config(DeltaRMax=0.5), [(5.0, 0.5, 0.49)])
    expected = c.Momentum + dressings[0].Momentum
    assert out.At(0).Momentum.E() == pytest.approx(expected.E(), rel=1e-12)


def test_dressing_just_outside_delta_r_excluded(load_delphes):
    c, _, out = run_with_refs(load_delphes, make_module_config(DeltaRMax=0.5), [(5.0, 0.5, 0.51)])
    assert out.At(0).Momentum.E() == c.Momentum.E()
    assert out.At(0).Momentum.Pt() == c.Momentum.Pt()


def test_dressing_pt_at_hard_boundary_excluded(load_delphes):
    c, _, out = run_with_refs(load_delphes, make_module_config(DeltaRMax=0.5), [(0.1, 0.5, 0.1)])
    assert out.At(0).Momentum.E() == c.Momentum.E()


def test_dressing_pt_just_above_boundary_included(load_delphes):
    c, dressings, out = run_with_refs(load_delphes, make_module_config(DeltaRMax=0.5), [(0.11, 0.5, 0.1)])
    expected = c.Momentum + dressings[0].Momentum
    assert out.At(0).Momentum.E() == pytest.approx(expected.E(), rel=1e-12)


def test_dressed_candidate_keeps_mother_reference(load_delphes):
    c, _, out = run_with_refs(load_delphes, make_module_config(DeltaRMax=0.5), [(5.0, 0.5, 0.1)])
    subs = out.At(0).GetCandidates()
    assert subs.GetEntries() == 1
    assert subs.At(0).GetUniqueID() == c.GetUniqueID()
