from conftest import assert_deterministic, build_config, candidate_snapshots


def make_module_config(**extra):
    return build_config(
        "BTagging",
        {"JetInputArray": "Delphes/inputJets", "BitNumber": 0, "EfficiencyFormula": [5, 0.8]},
        **extra,
    )


def test_b_tags_b_jet(run_tagging):
    jets = run_tagging(make_module_config(), [{"flavor": 5}])
    assert jets.GetEntries() == 1
    jet = jets.At(0)
    assert jet.BTag & 1 == 1


def test_empty_input(run_tagging):
    jets = run_tagging(make_module_config(), [])
    assert jets.GetEntries() == 0


def test_no_b_tag_light_jet(run_tagging):
    config = make_module_config(EfficiencyFormula=[5, 0.8, 0, 0.0])
    jets = run_tagging(config, [{"flavor": 1}])
    assert jets.GetEntries() == 1
    jet = jets.At(0)
    assert jet.BTag & 1 == 0


def test_b_tag_with_high_efficiency(run_tagging):
    jets = run_tagging(make_module_config(EfficiencyFormula=[5, 1.0]), [{"flavor": 5}])
    assert jets.At(0).BTag & 1 == 1


def test_bit_number_1(run_tagging):
    jets = run_tagging(make_module_config(BitNumber=1, EfficiencyFormula=[5, 1.0]), [{"flavor": 5}])
    assert jets.At(0).BTag & 1 == 0
    assert jets.At(0).BTag & 2 == 2


def test_c_tagging(run_tagging):
    jets = run_tagging(make_module_config(EfficiencyFormula=[4, 1.0]), [{"flavor": 4}])
    assert jets.At(0).BTag & 1 == 1


def test_deterministic_with_fixed_seed(run_tagging):
    config = make_module_config(EfficiencyFormula=[5, 0.5, 0, 0.0])
    jets_in = [{"flavor": 5} for _ in range(5)]
    assert_deterministic(
        lambda: run_tagging(config, jets_in),
        extract=lambda jets: candidate_snapshots(jets, ("BTag", "Flavor")),
    )
