from .conftest import (
    add_generator_particles,
    add_input_jets,
    add_input_partons,
    assert_deterministic,
    build_config,
    candidate_snapshots,
    run_seeds,
)


def make_module_config(**extra):
    return build_config(
        "TauTagging",
        {
            "ParticleInputArray": "Delphes/inputParticles",
            "PartonInputArray": "Delphes/inputPartons",
            "JetInputArray": "Delphes/inputJets",
            "BitNumber": 0,
            "DeltaR": 0.5,
            "TauPTMin": 0.0,
            "TauEtaMax": 10.0,
            "EfficiencyFormula": [15, 1.0],
        },
        **extra,
    )


def test_tau_tagging(run_tagging):
    jets = run_tagging(
        make_module_config(),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}],
        particles=[
            {"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1},
            {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1},
        ],
    )
    assert jets.GetEntries() == 1
    assert jets.At(0).TauTag & 1 == 1


def test_empty_input(run_tagging):
    jets = run_tagging(make_module_config(), jets=[], partons=[], particles=[])
    assert jets.GetEntries() == 0


def test_no_tau_tag_light_jet(run_tagging):
    jets = run_tagging(
        make_module_config(EfficiencyFormula=[15, 1.0, 0, 0.0]),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 0.5, "pid": 1, "status": 3, "d1": -1, "d2": -1}],
        particles=[],
    )
    assert jets.GetEntries() == 1
    assert jets.At(0).TauTag & 1 == 0


def test_bit_number_1(run_tagging):
    jets = run_tagging(
        make_module_config(BitNumber=1),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}],
        particles=[
            {"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1},
            {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1},
        ],
    )
    assert jets.GetEntries() == 1
    assert jets.At(0).TauTag & 1 == 0
    assert jets.At(0).TauTag & 2 == 2


def test_tau_pt_at_min_rejected(run_tagging):
    pions = [{"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1}]
    jets = run_tagging(
        make_module_config(TauPTMin=50.0),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.0, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}],
        particles=pions,
    )
    assert jets.At(0).TauTag & 1 == 0
    jets = run_tagging(
        make_module_config(TauPTMin=50.0),
        jets=[{"pt": 50.0, "eta": 0.5}],
        partons=[{"pt": 50.01, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}],
        particles=pions,
    )
    assert jets.At(0).TauTag & 1 == 1


def test_tau_eta_at_max_accepted(run_tagging):
    jets = run_tagging(
        make_module_config(TauEtaMax=2.5),
        jets=[{"pt": 50.0, "eta": 2.5}],
        partons=[{"pt": 50.0, "eta": 2.49, "pid": 15, "status": 3, "d1": 0, "d2": 1}],
        particles=[
            {"pt": 20.0, "eta": 2.4, "pid": 211, "status": 1},
            {"pt": 15.0, "eta": 2.6, "pid": -211, "status": 1},
        ],
    )
    assert jets.At(0).TauTag & 1 == 1
    jets = run_tagging(
        make_module_config(TauEtaMax=2.5),
        jets=[{"pt": 50.0, "eta": 2.5}],
        partons=[{"pt": 50.0, "eta": 2.51, "pid": 15, "status": 3, "d1": 0, "d2": 1}],
        particles=[
            {"pt": 20.0, "eta": 2.4, "pid": 211, "status": 1},
            {"pt": 15.0, "eta": 2.6, "pid": -211, "status": 1},
        ],
    )
    assert jets.At(0).TauTag & 1 == 0


def test_efficiency_zero_never_tags(load_delphes):
    config = make_module_config(EfficiencyFormula=[15, 0.0])

    def setup(module, factory):
        add_input_jets(module, factory, [{"pt": 50.0, "eta": 0.5}])
        add_input_partons(module, factory, [{"pt": 50.0, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}])
        add_generator_particles(
            module,
            factory,
            [{"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1}],
        )

    for jets in run_seeds(load_delphes, config, range(42, 62), setup, output="Delphes/inputJets"):
        assert jets.GetEntries() == 1
        assert jets.At(0).TauTag & 1 == 0


def test_efficiency_one_always_tags(load_delphes):
    config = make_module_config(EfficiencyFormula=[15, 1.0])

    def setup(module, factory):
        add_input_jets(module, factory, [{"pt": 50.0, "eta": 0.5}])
        add_input_partons(module, factory, [{"pt": 50.0, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}])
        add_generator_particles(
            module,
            factory,
            [{"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1}, {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1}],
        )

    for jets in run_seeds(load_delphes, config, range(42, 62), setup, output="Delphes/inputJets"):
        assert jets.At(0).TauTag & 1 == 1


def test_deterministic_with_fixed_seed(run_tagging):
    config = make_module_config(EfficiencyFormula=[15, 0.5])
    jets = [{"pt": 50.0, "eta": 0.5}]
    partons = [{"pt": 50.0, "eta": 0.5, "pid": 15, "status": 3, "d1": 0, "d2": 1}]
    particles = [
        {"pt": 20.0, "eta": 0.4, "pid": 211, "status": 1},
        {"pt": 15.0, "eta": 0.6, "pid": -211, "status": 1},
    ]
    assert_deterministic(
        lambda: run_tagging(config, jets=jets, partons=partons, particles=particles),
        extract=lambda out: candidate_snapshots(out, ("TauTag", "Charge")),
    )
