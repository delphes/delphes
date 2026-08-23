import cppyy
import pytest
import ROOT
from conftest import read_tree_branch

GET_CONFIGURED_CASES = [
    ("GetInt", "MyInt", 12345, -1),
    ("GetLong", "MyLong", -987654321, 0),
    ("GetDouble", "MyDouble", 0.25, -0.5),
    ("GetBool", "MyBool", False, True),
    ("GetString", "MyString", "hello_world", "sentinel"),
]


GET_DEFAULT_CASES = [
    ("GetInt", -11),
    ("GetLong", -22),
    ("GetDouble", -3.5),
    ("GetBool", True),
    ("GetString", "dflt"),
]


def module_scoped_config(**params):
    return {**{f"Delphes::{key}": value for key, value in params.items()}}


@pytest.mark.parametrize("method,default", GET_DEFAULT_CASES, ids=[case[0] for case in GET_DEFAULT_CASES])
def test_get_returns_default_when_key_absent(load_delphes, method, default):
    module, _, _ = load_delphes({})

    assert getattr(module, method)("NoSuchKey", default) == default

    assert getattr(module, method)("NoSuchKey", default, 0) == default


@pytest.mark.parametrize(
    "method,key,value,sentinel", GET_CONFIGURED_CASES, ids=[case[0] for case in GET_CONFIGURED_CASES]
)
def test_get_returns_configured_value(load_delphes, method, key, value, sentinel):
    module, _, _ = load_delphes(module_scoped_config(**{key: value}))

    assert getattr(module, method)(key, sentinel) == value


def test_get_param_missing_key_is_empty(load_delphes):
    module, _, _ = load_delphes({})
    param = module.GetParam("NoSuchKey")

    assert param.GetSize() == 0
    assert param.GetInt(7) == 7
    assert param.GetLong(-8) == -8
    assert param.GetDouble(8.5) == 8.5
    assert param.GetBool(False) is False
    assert param.GetString("dflt") == "dflt"


def test_get_param_list_surface(load_delphes):
    module, _, _ = load_delphes(module_scoped_config(MyList=[10, 20, 30]))
    param = module.GetParam("MyList")
    assert param.GetSize() == 3

    for i, expected in enumerate((10, 20, 30)):
        assert param[i].GetInt() == expected
        assert module.GetInt("MyList", -1, i) == expected

    assert module.GetParam("NoSuchKey")[0].GetSize() == 0


def test_get_string_numeric_string_stays_string(load_delphes):
    module, _, _ = load_delphes(module_scoped_config(NumStr="42"))
    assert module.GetString("NumStr", "dflt") == "42"
    assert module.GetInt("NumStr", -1) == 42


def test_export_array_is_not_idempotent_import_resolves_first(load_delphes):
    module, factory, _ = load_delphes({})

    first = module.ExportArray("dup")
    second = module.ExportArray("dup")
    assert cppyy.addressof(first) != cppyy.addressof(second)
    assert first.GetName() == "Delphes/dup"
    assert second.GetName() == "Delphes/dup"

    orphan = factory.NewCandidate()
    orphan.PID = 211
    second.Add(orphan)
    resolved = module.ImportArray("Delphes/dup")
    assert cppyy.addressof(resolved) == cppyy.addressof(first)
    assert resolved.GetEntries() == 0

    visible = factory.NewCandidate()
    visible.PID = 13
    first.Add(visible)
    assert module.ImportArray("Delphes/dup").GetEntries() == 1
    assert module.ImportArray("Delphes/dup").At(0).PID == 13


def test_new_branch_readable_after_write(load_delphes, tmp_path):
    def build(module, writer):
        branch = module.NewBranch("Extra", ROOT.Candidate.Class())
        tree = writer.GetTree()

        assert tree.GetBranch("Extra") is not None
        assert tree.GetBranch("Extra_size") is not None

        entry = branch.NewEntry()
        entry.PID = 211

    tree_reader, fin = read_tree_branch(load_delphes, tmp_path, build)

    branch = tree_reader.UseBranch("Extra")
    tree_reader.ReadEntry(0)
    assert branch.GetEntries() == 1
    assert branch.At(0).PID == 211

    fin.Close()


def test_add_info_readable_after_write(load_delphes, tmp_path):
    def build(module, writer):
        module.AddInfo("testInfo", 12.5)
        module.AddInfo("otherInfo", -1.25)
        user_info = writer.GetTree().GetUserInfo()
        assert user_info.GetSize() == 2
        assert user_info.At(0).GetName() == "testInfo"
        assert user_info.At(0).GetVal() == 12.5
        assert user_info.At(1).GetName() == "otherInfo"
        assert user_info.At(1).GetVal() == -1.25

    _, fin = read_tree_branch(load_delphes, tmp_path, build)

    user_info = fin.Get("Delphes").GetUserInfo()
    assert user_info.GetSize() == 2
    assert user_info.At(0).GetName() == "testInfo"
    assert user_info.At(0).GetVal() == 12.5
    assert user_info.At(1).GetName() == "otherInfo"
    assert user_info.At(1).GetVal() == -1.25

    fin.Close()


def test_add_info_and_new_branch_without_file_are_safe_noops(load_delphes):
    module, _, _ = load_delphes({})
    module.AddInfo("noFile", 1.0)
    branch = module.NewBranch("detached", ROOT.Candidate.Class())
    assert branch is not None

    module.GetTreeWriter().Fill()
    module.GetTreeWriter().Clear()
