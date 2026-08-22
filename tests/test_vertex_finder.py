from conftest import VERTEX_FINDER_TRACKS, build_config, run_vertex_finder_test


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
    tracks, vertices = run_vertex_finder_test(run_generic, make_module_config(), VERTEX_FINDER_TRACKS)
    assert vertices.GetEntries() >= 1
    assert vertices.At(0).ClusterNDF >= 4


def test_no_vertex_with_few_tracks(run_generic):
    tracks, vertices = run_vertex_finder_test(
        run_generic, make_module_config(), [(50.0, 0.5, 0.0, 0.001, 0), (40.0, 0.5, 0.001, 0.001, 0)]
    )
    assert vertices.GetEntries() == 0


def test_tracks_assigned_to_vertex(run_generic):
    tracks, vertices = run_vertex_finder_test(run_generic, make_module_config(), VERTEX_FINDER_TRACKS)
    assert tracks.GetEntries() == 4
    for i in range(tracks.GetEntries()):
        assert tracks.At(i).ClusterIndex >= 0
