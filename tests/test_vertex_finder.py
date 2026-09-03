from .conftest import VERTEX_FINDER_TRACKS, build_config, run_vertex_finder_test


def make_module_config(**extra):
    return build_config(
        "VertexFinder",
        {
            "InputArray": "Delphes/inputTracks",
            "OutputArray": "tracks",
            "VertexOutputArray": "vertices",
            "Sigma": 3.0,
            "MinPT": 0.0,
            "MaxEta": 10.0,
            "SeedMinPT": 0.0,
            "MinNDF": 4,
            "GrowSeeds": 1,
        },
        **extra,
    )


def test_finds_vertex_with_enough_tracks(run_generic):
    _, vertices = run_vertex_finder_test(run_generic, make_module_config(), VERTEX_FINDER_TRACKS)
    assert vertices.GetEntries() >= 1
    assert vertices.At(0).ClusterNDF >= 4


def test_no_vertex_with_few_tracks(run_generic):
    _, vertices = run_vertex_finder_test(
        run_generic, make_module_config(), [(50.0, 0.5, 0.0, 0.001, 0), (40.0, 0.5, 0.001, 0.001, 0)]
    )
    assert vertices.GetEntries() == 0


def test_tracks_assigned_to_vertex(run_generic):
    tracks, _ = run_vertex_finder_test(run_generic, make_module_config(), VERTEX_FINDER_TRACKS)
    assert tracks.GetEntries() == 4
    for i in range(tracks.GetEntries()):
        assert tracks.At(i).ClusterIndex >= 0


def test_vertex_position_accuracy(run_generic):
    _, vertices = run_vertex_finder_test(run_generic, make_module_config(), [(50.0, 0.5, 0.0, 0.001, 0)] * 4)
    assert vertices.GetEntries() == 1

    assert abs(vertices.At(0).Position.Z()) < 1e-6


def test_pileup_tracks_form_separate_cluster(run_generic):
    pv = [(50.0, 0.5, 0.0, 0.001, 0)] * 4
    pu = [(30.0, 0.5, 10.0, 0.001, 1)] * 4
    _, vertices = run_vertex_finder_test(run_generic, make_module_config(), pv + pu)
    assert vertices.GetEntries() == 2

    assert abs(vertices.At(0).Position.Z()) < 1e-6
    assert abs(vertices.At(1).Position.Z() - 10000.0) < 1e-3


def test_empty_input(run_generic):
    tracks, vertices = run_vertex_finder_test(run_generic, make_module_config(), [])
    assert tracks.GetEntries() == 0
    assert vertices.GetEntries() == 0
