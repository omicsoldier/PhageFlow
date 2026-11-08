"""Tests for configuration loading and database path precedence."""
from importlib.resources import files
from pathlib import Path

import pytest

from phageflow import __version__
from phageflow.utils.config import Config, load_config

BUNDLED = Path(str(files("phageflow.config").joinpath("default_config.yaml")))


@pytest.fixture(autouse=True)
def _no_env_db(monkeypatch):
    monkeypatch.delenv("PHAGEFLOW_DB", raising=False)


def test_dataclass_defaults():
    cfg = Config()
    assert cfg.version == __version__
    assert cfg.qc.average_qual == 20
    assert cfg.qc.length_required == 50
    assert cfg.genomad.composition == "auto"
    assert cfg.genomad.rescue_min_length_bp == 3000
    assert cfg.checkv.large_nd_rescue_bp == 20000
    assert "141" not in cfg.assembly.kmers
    assert cfg.threads >= 1


def test_bundled_config_matches_dataclass(tmp_path):
    """The shipped YAML and the dataclass defaults must not drift apart."""
    cfg = load_config(BUNDLED, workdir=tmp_path, is_bundled=True)
    ref = Config()
    for section in ("qc", "assembly", "host_removal", "genomad", "checkv", "annotate"):
        assert getattr(cfg, section) == getattr(ref, section), section


def test_empty_config_file(tmp_path):
    p = tmp_path / "config" / "config.yaml"
    p.parent.mkdir()
    p.write_text("")
    cfg = load_config(p)
    assert cfg.workdir == tmp_path
    assert cfg.databases.genomad == tmp_path / "databases" / "genomad_db"


def test_null_sections_are_tolerated(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("qc:\nhost_removal:\n  always_include_accessions:\ndatabases:\n")
    cfg = load_config(p, workdir=tmp_path)
    assert cfg.host_removal.always_include_accessions == []


def test_relative_and_absolute_db_paths(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(f"databases:\n  kraken2: k2\n  phold: {tmp_path / 'abs_phold'}\n")
    cfg = load_config(p, workdir=tmp_path)
    assert cfg.databases.kraken2 == tmp_path / "k2"
    assert cfg.databases.phold == tmp_path / "abs_phold"


def _make_env_db(root: Path) -> Path:
    (root / "checkv_db" / "checkv-db-v1.6").mkdir(parents=True)
    for d in ("genomad_db", "pharokka_db", "phold_db", "k2_db"):
        (root / d).mkdir()
    return root


def test_env_db_overrides_bundled_defaults(tmp_path, monkeypatch):
    db = _make_env_db(tmp_path / "dbs")
    monkeypatch.setenv("PHAGEFLOW_DB", str(db))
    cfg = load_config(BUNDLED, workdir=tmp_path, is_bundled=True)
    assert cfg.databases.genomad == db / "genomad_db"
    assert cfg.databases.kraken2 == db / "k2_db"
    # newest versioned CheckV directory is picked up
    assert cfg.databases.checkv == db / "checkv_db" / "checkv-db-v1.6"


def test_explicit_config_beats_env_db(tmp_path, monkeypatch):
    db = _make_env_db(tmp_path / "dbs")
    monkeypatch.setenv("PHAGEFLOW_DB", str(db))
    p = tmp_path / "c.yaml"
    p.write_text("databases:\n  genomad: /custom/genomad\n")
    cfg = load_config(p, workdir=tmp_path)
    assert cfg.databases.genomad == Path("/custom/genomad")
    # keys missing from the file fall back to PHAGEFLOW_DB
    assert cfg.databases.phold == db / "phold_db"


def test_section_overrides(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(
        "threads: 3\n"
        "genomad:\n  min_virus_hallmarks: 2\n  min_score: 0.8\n"
        "annotate:\n  phold_gpu: false\n"
    )
    cfg = load_config(p, workdir=tmp_path)
    assert cfg.threads == 3
    assert cfg.genomad.min_hallmarks == 2
    assert cfg.genomad.min_score == 0.8
    assert cfg.annotate.phold_gpu is False


def test_results_paths(tmp_path):
    cfg = Config(workdir=tmp_path)
    cfg.set_output_dir(tmp_path / "out" / "results")
    assert cfg.results("01_qc") == tmp_path / "out" / "results" / "01_qc"
    assert cfg.reports_dir == tmp_path / "reports"
