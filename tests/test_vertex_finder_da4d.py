from conftest import VERTEX_FINDER_TRACKS, build_config, make_vertex_finder_track, run_vertex_finder_test


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
    _, vertices = run_vertex_finder_test(run_generic, make_module_config(), VERTEX_FINDER_TRACKS)
    assert vertices.GetEntries() >= 1


def test_single_track_vertex(run_generic):
    _, vertices = run_vertex_finder_test(run_generic, make_module_config(), [(50.0, 0.5, 0.0, 0.001, 0)])
    assert vertices.GetEntries() == 1
    assert vertices.At(0).ClusterNDF == 1


def test_empty_input(run_generic):
    tracks, vertices = run_vertex_finder_test(run_generic, make_module_config(), [])
    assert tracks.GetEntries() == 0
    assert vertices.GetEntries() == 0


def test_repeated_process(load_delphes):
    config = make_module_config()
    module, factory, _ = load_delphes(config)
    input_array = module.ExportArray("inputTracks")
    for pt, eta, dz, error_dz, is_pu in VERTEX_FINDER_TRACKS:
        input_array.Add(make_vertex_finder_track(factory, pt, eta, dz, error_dz, is_pu))
    module.Init()
    module.Process()
    vertices = module.ImportArray("TestModule/vertices")
    first = vertices.GetEntries()
    assert first >= 1
    module.Process()
    assert vertices.GetEntries() >= 2
