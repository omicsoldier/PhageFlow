"""Runner orchestration tests with all modules stubbed out."""
import sys
import types
from pathlib import Path

import pytest

from phageflow.modules import runner
from phageflow.utils import status
from phageflow.utils.config import Config


@pytest.fixture
def cfg(tmp_path: Path) -> Config:
    c = Config(workdir=tmp_path)
    c.set_output_dir(tmp_path / "results")
    c.set_reports_dir(tmp_path / "reports")
    return c


def _stub(monkeypatch, name: str, fn) -> None:
    mod = types.ModuleType(f"phageflow.modules.{name}")
    mod.run = fn
    monkeypatch.setitem(sys.modules, f"phageflow.modules.{name}", mod)


def _install_stubs(monkeypatch, tmp_path: Path, fail_on: str | None = None):
    calls: list[tuple[str, str]] = []

    def make(module, ret):
        def fn(cfg, sample_id, **kw):
            calls.append((module, sample_id))
            if module == fail_on:
                raise RuntimeError(f"{module} boom")
            return ret(sample_id, kw)
        return fn

    def touch(p: Path) -> Path:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("x")
        return p

    def ann_ready(sid, kw):
        d = tmp_path / "ann" / sid
        touch(d / "Phage_candidate_001.fasta")
        return d, tmp_path / "drafts"

    _stub(monkeypatch, "qc", make("qc", lambda s, kw: (touch(tmp_path / f"{s}_qc1"), touch(tmp_path / f"{s}_qc2"))))
    _stub(monkeypatch, "host_removal", make(
        "host-removal", lambda s, kw: (touch(tmp_path / f"{s}_h1"), touch(tmp_path / f"{s}_h2"), tmp_path / "none")))
    _stub(monkeypatch, "assembly", make("assembly", lambda s, kw: touch(tmp_path / f"{s}_contigs.fasta")))
    _stub(monkeypatch, "coverage", make("coverage", lambda s, kw: None))
    _stub(monkeypatch, "viral_id", make("viral-id", lambda s, kw: (touch(tmp_path / f"{s}_virus.fna"), None)))
    _stub(monkeypatch, "quality", make("quality", ann_ready))
    _stub(monkeypatch, "annotate", make("annotate", lambda s, kw: None))
    _stub(monkeypatch, "resistance", make("resistance", lambda s, kw: None))
    return calls


def test_full_run_records_status(cfg, tmp_path, monkeypatch):
    calls = _install_stubs(monkeypatch, tmp_path)
    pairs = [(tmp_path / "A_R1.fastq.gz", tmp_path / "A_R2.fastq.gz")]
    res = runner.run(cfg, pairs)
    assert res == {"A": "ok"}
    assert [m for m, _ in calls] == runner._MODULE_ORDER
    rows = {r["module"]: r["status"] for r in status.get_all(cfg.reports_dir)}
    assert rows == {m: "completed" for m in runner._MODULE_ORDER}


def test_failure_isolated_per_sample(cfg, tmp_path, monkeypatch):
    _install_stubs(monkeypatch, tmp_path, fail_on="assembly")
    pairs = [
        (tmp_path / "A_R1.fastq.gz", tmp_path / "A_R2.fastq.gz"),
        (tmp_path / "B_R1.fastq.gz", tmp_path / "B_R2.fastq.gz"),
    ]
    res = runner.run(cfg, pairs)
    assert set(res) == {"A", "B"}
    assert all(v.startswith("failed") for v in res.values())
    rows = [(r["sample_id"], r["module"], r["status"]) for r in status.get_all(cfg.reports_dir)]
    assert ("A", "assembly", "failed") in rows and ("B", "assembly", "failed") in rows


def test_system_exit_does_not_abort_other_samples(cfg, tmp_path, monkeypatch):
    _install_stubs(monkeypatch, tmp_path)

    def exits(cfg, sample_id, **kw):
        if sample_id == "A":
            sys.exit(1)  # e.g. require_tools() with a missing binary
        return tmp_path / "x1", tmp_path / "x2"

    _stub(monkeypatch, "qc", exits)
    pairs = [
        (tmp_path / "A_R1.fastq.gz", tmp_path / "A_R2.fastq.gz"),
        (tmp_path / "B_R1.fastq.gz", tmp_path / "B_R2.fastq.gz"),
    ]
    res = runner.run(cfg, pairs)
    assert res["A"].startswith("failed")
    assert "B" in res


def test_coverage_failure_is_non_fatal(cfg, tmp_path, monkeypatch):
    _install_stubs(monkeypatch, tmp_path, fail_on="coverage")
    res = runner.run(cfg, [(tmp_path / "A_R1.fastq.gz", tmp_path / "A_R2.fastq.gz")])
    assert res == {"A": "ok"}
    rows = {r["module"]: r["status"] for r in status.get_all(cfg.reports_dir)}
    assert rows["coverage"] == "failed"
    assert rows["resistance"] == "completed"


def test_duplicate_sample_ids_rejected(cfg, tmp_path):
    pairs = [
        (tmp_path / "a" / "S_R1.fastq.gz", tmp_path / "a" / "S_R2.fastq.gz"),
        (tmp_path / "b" / "S_R1.fastq.gz", tmp_path / "b" / "S_R2.fastq.gz"),
    ]
    with pytest.raises(ValueError, match="same sample ID"):
        runner.run(cfg, pairs)


def test_resume_missing_input_fails_clearly(cfg, tmp_path, monkeypatch):
    _install_stubs(monkeypatch, tmp_path)
    res = runner.run(cfg, [(tmp_path / "A_R1.fastq.gz", tmp_path / "A_R2.fastq.gz")],
                     from_module="viral-id")
    assert res["A"].startswith("failed") and "Resume input" in res["A"]
