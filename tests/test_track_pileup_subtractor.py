import pytest
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


def test_empty_input(run_generic):
    def setup(module, factory):
        module.ExportArray("inputVertices")
        module.ExportArray("inputTracks")

    output = run_generic(make_module_config(), setup=setup, outputs=("TestModule/outputTracks",))
    assert output.GetEntries() == 0


def run_multi_test(run_generic, config, vertex_specs, track_specs):
    def setup(module, factory):
        vertex_array = module.ExportArray("inputVertices")
        for z, is_pu in vertex_specs:
            v = make_vertex(factory, z=z, is_pu=is_pu)
            vertex_array.Add(v)

        track_array = module.ExportArray("inputTracks")
        for pt, z, is_pu, charge in track_specs:
            t = make_candidate(factory, pt, 0.5, pid=211, charge=charge)
            t.Position.SetXYZT(0.0, 0.0, z, 0.0)
            t.IsPU = is_pu
            mother = make_candidate(factory, pt, 0.5, pid=211, charge=charge)
            mother.Position.SetXYZT(0.0, 0.0, z, 0.0)
            t.AddCandidate(mother)
            track_array.Add(t)

    return run_generic(config, setup=setup, outputs=("TestModule/outputTracks",))


def test_dz_resolution_boundary(run_generic):
    kept = run_multi_test(run_generic, make_module_config(), [(0.0, 0)], [(50.0, 1.0, 1, 1)])
    assert kept.GetEntries() == 1
    assert kept.At(0).IsRecoPU == 0
    assert kept.At(0).IsPU == 1
    dropped = run_multi_test(run_generic, make_module_config(), [(0.0, 0)], [(50.0, 1.0001, 1, 1)])
    assert dropped.GetEntries() == 0


def test_resolution_scales_with_formula(run_generic):
    tight = run_multi_test(run_generic, make_module_config(), [(0.0, 0)], [(50.0, 5.0, 1, 1)])
    assert tight.GetEntries() == 0
    wide = run_multi_test(run_generic, make_module_config(ZVertexResolution=0.010), [(0.0, 0)], [(50.0, 5.0, 1, 1)])
    assert wide.GetEntries() == 1
    assert wide.At(0).IsRecoPU == 0


def test_ptmin_strict_boundary(run_generic):
    dropped = run_multi_test(run_generic, make_module_config(PTMin=50.0), [(0.0, 0)], [(50.0, 0.0, 0, 1)])
    assert dropped.GetEntries() == 0
    kept = run_multi_test(run_generic, make_module_config(PTMin=50.0), [(0.0, 0)], [(50.01, 0.0, 0, 1)])
    assert kept.GetEntries() == 1
    assert kept.At(0).Momentum.Pt() == pytest.approx(50.01, rel=1e-9)


def test_neutral_pu_track_is_kept(run_generic):
    output = run_multi_test(run_generic, make_module_config(), [(0.0, 0)], [(50.0, 100.0, 1, 0)])
    assert output.GetEntries() == 1
    assert output.At(0).IsRecoPU == 0
    assert output.At(0).IsPU == 1


def test_last_non_pu_vertex_defines_pv(run_generic):
    output = run_multi_test(run_generic, make_module_config(), [(10.0, 0), (20.0, 0)], [(50.0, 11.0, 1, 1)])
    assert output.GetEntries() == 0

    output = run_multi_test(run_generic, make_module_config(), [(10.0, 0), (15.0, 1), (20.0, 0)], [(50.0, 11.0, 1, 1)])
    assert output.GetEntries() == 0

    output = run_multi_test(run_generic, make_module_config(), [(20.0, 0), (10.0, 0)], [(50.0, 11.0, 1, 1)])
    assert output.GetEntries() == 1
    assert output.At(0).IsRecoPU == 0
