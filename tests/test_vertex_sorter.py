import pytest
from .conftest import build_config, make_candidate, make_vertex


def make_module_config(**extra):
    return build_config(
        "VertexSorter",
        {
            "InputArray": "Delphes/inputVertices",
            "TrackInputArray": "Delphes/inputTracks",
            "OutputArray": "clusters",
            "Method": "GenBest",
        },
        **extra,
    )


def run_sorter_test(run_generic, config, vertexs, tracks):
    def setup(module, factory):
        vertex_array = module.ExportArray("inputVertices")
        for z, ndf, sumpt2, cluster_index in vertexs:
            v = make_vertex(factory, z=z * 1000.0)
            v.ClusterNDF = ndf
            v.SumPT2 = sumpt2
            v.ClusterIndex = cluster_index
            vertex_array.Add(v)

        track_array = module.ExportArray("inputTracks")
        for pt, eta, cluster_index, is_pu in tracks:
            t = make_candidate(factory, pt, eta, charge=1, pid=211)
            t.ClusterIndex = cluster_index
            t.IsPU = is_pu
            track_array.Add(t)

    return run_generic(config, setup=setup, outputs=("TestModule/clusters",))


def test_sorts_by_gen_best(run_generic):
    output = run_sorter_test(
        run_generic,
        make_module_config(),
        vertexs=[(0.0, 4, 100.0, 0), (1.0, 4, 200.0, 1)],
        tracks=[(10.0, 0.5, 0, 0), (20.0, 0.5, 1, 0)],
    )
    assert output.GetEntries() == 2
    assert output.At(0).GenSumPT2 == pytest.approx(400.0, rel=1e-3)
    assert output.At(1).GenSumPT2 == pytest.approx(100.0, rel=1e-3)


def test_single_vertex(run_generic):
    output = run_sorter_test(
        run_generic,
        make_module_config(),
        vertexs=[(0.0, 4, 150.0, 0)],
        tracks=[(10.0, 0.5, 0, 0), (20.0, 0.5, 0, 0)],
    )
    assert output.GetEntries() == 1
    assert output.At(0).GenSumPT2 == pytest.approx(500.0, rel=1e-3)


def test_empty_input(run_generic):
    output = run_sorter_test(run_generic, make_module_config(), vertexs=[], tracks=[])
    assert output.GetEntries() == 0
