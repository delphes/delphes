from collections import Counter

import pytest
from conftest import assert_deterministic, build_config, candidate_snapshots, make_candidate, run_repeated


def make_module_config(**extra):
    return build_config("IdentificationMap", {"EfficiencyFormula": [13, 11, 1.0]}, **extra)


def test_remapping(run_module):
    output = run_module(make_module_config(), [{"pid": 13, "charge": 1}])
    assert output.GetEntries() == 1
    assert output.At(0).PID == 11


def test_no_match_passes_through(run_module):
    output = run_module(make_module_config(), [{"pid": 211, "charge": 1}])
    assert output.GetEntries() == 1
    assert output.At(0).PID == 211


def test_zero_efficiency_rejects(run_module):
    output = run_module(make_module_config(EfficiencyFormula=[13, 11, 0.0]), [{"pid": 13, "charge": 1}])
    assert output.GetEntries() == 0


def test_empty_input(run_module):
    output = run_module(make_module_config(), [])
    assert output.GetEntries() == 0


def test_full_mapping_table(run_module):
    config = make_module_config(EfficiencyFormula=[13, 11, 1.0, 211, 2212, 1.0, 11, 13, 1.0])
    output = run_module(config, [{"pid": 13, "charge": 1}, {"pid": 211, "charge": 1}, {"pid": 11, "charge": -1}])
    assert output.GetEntries() == 3
    assert [output.At(i).PID for i in range(3)] == [11, 2212, -13]


def test_neutral_particle_remapped_to_pid_zero(run_module):
    config = make_module_config(EfficiencyFormula=[22, 23, 1.0])
    output = run_module(config, [{"pid": 22, "charge": 0}])
    assert output.GetEntries() == 1
    assert output.At(0).PID == 0


def test_negative_input_pid_falls_back_to_positive_entry(run_module):
    output = run_module(make_module_config(), [{"pid": -13, "charge": -1}])
    assert output.GetEntries() == 1
    assert output.At(0).PID == -11


def test_remapped_candidate_keeps_momentum(run_module):
    output = run_module(make_module_config(), [{"pt": 37.0, "eta": 1.25, "phi": 0.7, "pid": 13, "charge": 1}])
    assert output.GetEntries() == 1
    c = output.At(0)
    assert c.Momentum.Pt() == pytest.approx(37.0, rel=1e-9)
    assert c.Momentum.Eta() == pytest.approx(1.25, rel=1e-9)
    assert c.Momentum.Phi() == pytest.approx(0.7, rel=1e-9)


def test_output_pid_zero_keeps_original_pid(run_module):
    config = make_module_config(EfficiencyFormula=[22, 0, 1.0])
    output = run_module(config, [{"pid": 22, "charge": 0}])
    assert output.GetEntries() == 1
    assert output.At(0).PID == 22


def test_default_entry_overridable(run_module):
    config = make_module_config(EfficiencyFormula=[0, 11, 1.0])
    output = run_module(config, [{"pid": 211, "charge": 1}])
    assert output.GetEntries() == 1
    assert output.At(0).PID == 11


def test_default_entry_zero_efficiency_drops_unmatched(run_module):
    config = make_module_config(EfficiencyFormula=[0, 0, 0.0])
    output = run_module(config, [{"pid": 211, "charge": 1}])
    assert output.GetEntries() == 0


def test_multimap_marginal_probabilities(load_delphes):
    config = make_module_config(EfficiencyFormula=[13, 11, 0.5, 13, 23, 0.5])

    def setup(module, factory):
        array = module.ExportArray("inputParticles")
        array.Add(make_candidate(factory, 50.0, 0.5, pid=13, charge=1))

    results = run_repeated(load_delphes, config, 400, setup)
    counts = Counter()
    for result in results:
        assert result.GetEntries() == 1
        counts[result.At(0).PID] += 1
    assert set(counts) == {11, 23}
    assert counts[11] == pytest.approx(200, abs=30)
    assert counts[23] == pytest.approx(200, abs=30)


def test_multimap_total_efficiency_below_one_drops(load_delphes):
    config = make_module_config(EfficiencyFormula=[13, 11, 0.3, 13, 23, 0.3])

    def setup(module, factory):
        array = module.ExportArray("inputParticles")
        array.Add(make_candidate(factory, 50.0, 0.5, pid=13, charge=1))

    results = run_repeated(load_delphes, config, 400, setup)
    kept = sum(1 for r in results if r.GetEntries() == 1)
    assert kept == pytest.approx(240, abs=40)


def test_deterministic_with_fixed_seed(run_generic):
    config = make_module_config(EfficiencyFormula=[13, 11, 0.5, 13, 23, 0.5])

    def setup(module, factory):
        array = module.ExportArray("inputParticles")
        for i in range(2):
            array.Add(make_candidate(factory, 50.0, 0.5 + 0.1 * i, pid=13, charge=1))

    assert_deterministic(
        lambda: run_generic(config, setup=setup),
        extract=lambda out: candidate_snapshots(out, ("PID",)),
    )
