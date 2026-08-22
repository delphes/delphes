from conftest import build_config, make_candidate, make_vertex


def make_module_config(**extra):
    return build_config(
        "TrackPileUpSubtractor",
        {
            "VertexInputArray": "Delphes/inputVertices",
            "ZVertexResolution": 0.001,
            "PTMin": 0.0,
            "InputArray": ["Delphes/inputTracks", "outputTracks"],
        },
        **extra,
    )


def run_subtractor_test(run_generic, config, vertex_z=0.0, vertex_is_pu=0, track_z=0.0, track_is_pu=0, track_charge=1):
    def setup(module, factory):
        vertex_array = module.ExportArray("inputVertices")
        v = make_vertex(factory, z=vertex_z * 1000.0, is_pu=vertex_is_pu)
        vertex_array.Add(v)

        track_array = module.ExportArray("inputTracks")
        t = make_candidate(factory, 50.0, 0.5, pid=211, charge=track_charge)
        t.Position.SetXYZT(0.0, 0.0, track_z * 1000.0, 0.0)
        t.IsPU = track_is_pu

        mother = make_candidate(factory, 50.0, 0.5, pid=211, charge=track_charge)
        mother.Position.SetXYZT(0.0, 0.0, track_z * 1000.0, 0.0)
        t.AddCandidate(mother)

        track_array.Add(t)

    return run_generic(config, setup=setup, outputs=("TestModule/outputTracks",))


def test_keeps_signal_track(run_generic):
    output = run_subtractor_test(run_generic, make_module_config(), vertex_z=0.0, track_z=0.0, track_is_pu=0)
    assert output.GetEntries() == 1
    assert output.At(0).IsRecoPU == 0


def test_subtracts_pileup_track(run_generic):
    output = run_subtractor_test(run_generic, make_module_config(), vertex_z=0.0, track_z=1.0, track_is_pu=1)
    assert output.GetEntries() == 0
