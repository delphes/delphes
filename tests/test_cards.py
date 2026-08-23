from pathlib import Path

import pytest

tests_dir = Path(__file__).parent
cards_dir = tests_dir.parent / "cards"

keyword = "ExecutionPath"

files = sorted(p for p in cards_dir.rglob("*.tcl") if p.is_file() and keyword in p.read_text(encoding="utf-8"))


@pytest.mark.parametrize("card_file", files, ids=lambda p: str(p.relative_to(cards_dir)))
def test_cards(load_delphes, card_file):
    module, _, _ = load_delphes(str(cards_dir / card_file))
    for name in ("allParticles", "stableParticles", "partons"):
        module.ExportArray(name)
    module.Init()
