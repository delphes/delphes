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
