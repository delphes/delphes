from conftest import VERTEX_FINDER_TRACKS, build_config, run_vertex_finder_test


def make_module_config(**extra):
    return build_config(
        "VertexFinderDA4D",
        {
            "InputArray": "Delphes/inputTracks",
            "OutputArray": "tracks",
            "VertexOutputArray": "vertices",
            "MinPT": 0.0,
            "VertexSpaceSize": 0.5,
            "VertexTimeSize": 10.0,
            "UseTc": 0,
            "BetaMax": 0.1,
            "BetaStop": 1.0,
            "CoolingFactor": 0.8,
            "MaxIterations": 100,
            "DzCutOff": 40.0,
            "D0CutOff": 30.0,
            "DtCutOff": 100.0,
        },
        **extra,
    )


def test_finds_vertex(run_generic):
    tracks, vertices = run_vertex_finder_test(run_generic, make_module_config(), VERTEX_FINDER_TRACKS)
    assert vertices.GetEntries() >= 1


def test_single_track_vertex(run_generic):
    tracks, vertices = run_vertex_finder_test(run_generic, make_module_config(), [(50.0, 0.5, 0.0, 0.001, 0)])
    assert vertices.GetEntries() == 1
    assert vertices.At(0).ClusterNDF == 1
