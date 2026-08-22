from conftest import build_config, make_candidate


def make_module_config(**extra):
    return build_config(
        "PhotonConversions",
        {
            "Radius": 1.0,
            "HalfLength": 3.0,
            "EtaMin": 2.0,
            "EtaMax": 5.0,
            "Step": 0.1,
            "ConversionMap": 0.0,
        },
        **extra,
    )


def run_conversion_test(run_generic, config, pid=22, pt=50.0, eta=3.0):
    def setup(module, factory):
        input_array = module.ExportArray("inputParticles")
        c = make_candidate(factory, pt, eta, pid=pid)
        input_array.Add(c)

    return run_generic(config, setup=setup, outputs=("TestModule/outputParticles",))


def test_non_photon_passes_through(run_generic):
    output = run_conversion_test(run_generic, make_module_config(), pid=211, pt=50.0, eta=3.0)
    assert output.GetEntries() == 1
    assert output.At(0).PID == 211


def test_photon_with_zero_conversion(run_generic):
    output = run_conversion_test(run_generic, make_module_config(), pid=22, pt=50.0, eta=3.0)
    assert output.GetEntries() == 1
    assert output.At(0).PID == 22
