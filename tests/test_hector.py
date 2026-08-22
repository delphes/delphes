from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "Hector",
        {
            "OutputArray": "hits",
            "Direction": 1,
            "BeamLineLength": 430.0,
            "Distance": 420.0,
            "OffsetX": 0.0,
            "OffsetS": 120.0,
            "SigmaE": 0.0,
            "SigmaX": 0.0,
            "SigmaY": 0.0,
            "SigmaT": 0.0,
            "EtaMin": 5.0,
            "BeamLineFile": "cards/LHCB1IR5_5TeV.tfs",
            "IPName": "IP5",
        },
        **extra,
    )


def run_hector_test(run_generic, config, eta):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, 50.0, eta, pid=211, charge=1)
        input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/hits",))


def test_forward_particle_propagates(run_generic):
    hits = run_hector_test(run_generic, make_module_config(), eta=6.0)
    assert hits.GetEntries() == 1
    assert hits.At(0).Momentum.E() > 0


def test_low_eta_particle_rejected(run_generic):
    hits = run_hector_test(run_generic, make_module_config(), eta=0.5)
    assert hits.GetEntries() == 0
