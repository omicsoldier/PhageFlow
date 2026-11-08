"""Tests for output parsers in the pipeline modules (no external tools needed)."""
import json
from pathlib import Path

from phageflow.modules import coverage, host_removal, qc, viral_id
from phageflow.utils import status
from phageflow.utils.tools import fasta_stats, q


# ── coverage ──────────────────────────────────────────────────────────────────

def test_coverm_header_trimmed_mean_not_confused_with_mean(tmp_path: Path):
    tsv = tmp_path / "cov.tsv"
    prefix = "s_contigs.fasta/s_R1.fastq.gz"
    tsv.write_text(
        f"Contig\t{prefix} Mean\t{prefix} Trimmed Mean\t{prefix} Covered Bases\t{prefix} Variance\n"
        "c1\t100.0\t90.0\t5000\t400.0\n"
    )
    (m,) = coverage._parse_coverage(tsv)
    assert m["mean_cov"] == 100.0
    assert m["trimmed_mean"] == 90.0
    assert m["covered_bases"] == 5000
    assert abs(m["cv"] - 0.2) < 1e-9


# ── qc ────────────────────────────────────────────────────────────────────────

def test_fastp_json_and_summary_tsv(tmp_path: Path):
    j = tmp_path / "fastp.json"
    j.write_text(json.dumps({
        "summary": {
            "before_filtering": {"total_reads": 1_000_000},
            "after_filtering": {"total_reads": 900_000, "total_bases": 5_000_000,
                                "gc_content": 0.5, "q20_rate": 0.98, "q30_rate": 0.9},
        },
        "duplication": {"rate": 0.1},
    }))
    m = qc._parse_fastp_json(j)
    assert m["pct_pass"] == "90.0%"
    assert m["est_coverage_typical"] == "100x"
    out = tmp_path / "qc_summary.tsv"
    qc._save_tsv("s1", m, out)
    header, row = out.read_text().splitlines()
    rec = dict(zip(header.split("\t"), row.split("\t")))
    assert rec["est_coverage_50kb"] == "100x"


def test_read_length_detection(tmp_path: Path):
    fq = tmp_path / "r.fastq"
    fq.write_text("".join(f"@r{i}\n{'A' * 250}\n+\n{'I' * 250}\n" for i in range(10)))
    assert qc._detect_read_length(fq) == 250
    assert qc._read_length_params(250)["length_required"] == 75
    assert qc._read_length_params(150) == {}


# ── host removal (Kraken2 reports) ────────────────────────────────────────────

K2_REPORT = "\n".join([
    " 10.00\t100\t100\tU\t0\tunclassified",
    " 90.00\t900\t0\tR\t1\troot",
    " 60.00\t600\t0\tD\t2\t  Bacteria",
    " 60.00\t600\t0\tP\t1224\t    Pseudomonadota",
    " 50.00\t500\t0\tG\t561\t      Escherichia",
    " 50.00\t500\t500\tS\t562\t        Escherichia coli",
    "  5.00\t50\t50\tS\t573\t        Klebsiella pneumoniae",
    " 30.00\t300\t0\tD\t10239\t  Viruses",
    " 30.00\t300\t300\tS\t10710\t    Escherichia phage lambda",
]) + "\n"


def test_contaminant_taxids_are_species_only(tmp_path: Path):
    rep = tmp_path / "k2.report"
    rep.write_text(K2_REPORT)
    taxids = host_removal._parse_contaminants_from_report(rep, min_pct=1.0)
    # phylum (1224) and genus (561) must not be downloaded as references
    assert taxids == [562, 573]


def test_viral_taxids(tmp_path: Path):
    rep = tmp_path / "k2.report"
    rep.write_text(K2_REPORT)
    assert host_removal._parse_viral_taxids(rep) == {10239, 10710}


def test_level_a_genus_detection(tmp_path: Path, monkeypatch):
    from phageflow.utils.config import Config

    k2_db = tmp_path / "k2"; k2_db.mkdir()
    cfg = Config(workdir=tmp_path)
    cfg.databases.kraken2 = k2_db
    r1 = tmp_path / "r1.fq.gz"; r1.write_bytes(b"x")
    r2 = tmp_path / "r2.fq.gz"; r2.write_bytes(b"x")
    rpt = tmp_path / "rpt"; rpt.mkdir()

    def fake_run(cmd, **kw):
        if cmd[0] == "bash":  # seqtk subsample
            (rpt / "s_la_sub_R1.fq").write_text("@a\nA\n+\nI\n")
            (rpt / "s_la_sub_R2.fq").write_text("@a\nA\n+\nI\n")
        else:  # kraken2
            (rpt / "s_levelA_contamination.report").write_text(K2_REPORT)

    monkeypatch.setattr(host_removal.subprocess, "run", fake_run)
    warnings: list[str] = []
    host_removal._level_a_check(r1, r2, rpt, "s", cfg, warnings)
    assert len(warnings) == 1
    assert "Escherichia (50.0%)" in warnings[0]


# ── viral id ──────────────────────────────────────────────────────────────────

def test_virus_summary_tiers(tmp_path: Path):
    tsv = tmp_path / "x_virus_summary.tsv"
    tsv.write_text(
        "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\t"
        "fdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
        "SPAdes|NODE_1\t40000\tDTR\tNA\t50\t11\t0.95\tNA\t5\t10\t"
        "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;Drexlerviridae\n"
        "MEGAHIT|k141_2\t5000\tNo terminal repeats\tNA\t5\t15\t0.5\tNA\t0\t1\t\n"
        "MEGAHIT|k141_3\t1000\tNo terminal repeats\tNA\t1\t11\t0.5\tNA\t0\t1\t\n"
    )
    hits = viral_id._parse_virus_summary(tsv, 0.7, 0.4, 3000)
    assert hits["SPAdes|NODE_1"]["status"] == "main"
    assert hits["SPAdes|NODE_1"]["family"] == "Drexlerviridae"
    assert hits["MEGAHIT|k141_2"]["status"] == "rescued"
    assert hits["MEGAHIT|k141_2"]["genetic_code"] == 15
    assert "MEGAHIT|k141_3" not in hits  # too short to rescue


def test_genetic_codes_from_genes_handles_pipes(tmp_path: Path):
    tsv = tmp_path / "x_virus_genes.tsv"
    tsv.write_text(
        "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\n"
        "SPAdes|NODE_1_length_100_1\t1\t90\t90\t1\t0.5\t15\n"
        "SPAdes|NODE_1_length_100_2\t100\t190\t90\t1\t0.5\t15\n"
        "other_1\t1\t90\t90\t1\t0.5\t11\n"
    )
    codes = viral_id._parse_genetic_codes_from_genes(tsv, {"SPAdes|NODE_1_length_100"})
    assert codes == {"SPAdes|NODE_1_length_100": 15}


def test_genomad_min_score_uses_rescue_threshold(tmp_path: Path, monkeypatch):
    from phageflow.utils.config import Config

    captured = {}
    monkeypatch.setattr(viral_id, "run_silent", lambda cmd, **kw: captured.setdefault("cmd", cmd))
    cfg = Config(workdir=tmp_path)
    contigs = tmp_path / "s_contigs.fasta"; contigs.write_text(">a\nACGT\n")
    viral_id._run_genomad(contigs, tmp_path / "g", tmp_path / "db", cfg, tmp_path / "log", True)
    cmd = captured["cmd"]
    assert cmd[cmd.index("--min-score") + 1] == str(cfg.genomad.rescue_min_score)


def test_parse_taxonomy():
    t = viral_id._parse_taxonomy("Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes")
    assert t == {"family": "Unknown", "finest": "Caudoviricetes", "naming_level": "Caudoviricetes"}
    assert viral_id._parse_taxonomy("")["naming_level"] == "Phage"


# ── utils ─────────────────────────────────────────────────────────────────────

def test_fasta_stats(tmp_path: Path):
    fa = tmp_path / "a.fa"
    fa.write_text(">a\nAAAA\nAA\n>b\nAAAA\n>c\nAA\n")
    assert fasta_stats(fa) == {"n": 3, "total_bp": 12, "largest_bp": 6, "n50": 6}


def test_shell_quote():
    assert q(Path("/a b/c")) == "'/a b/c'"


def test_status_roundtrip(tmp_path: Path):
    status.mark(tmp_path, "s1", "qc", "running")
    status.mark(tmp_path, "s1", "qc", "completed", notes="ok")
    status.mark(tmp_path, "s2", "qc", "failed")
    assert status.is_completed(tmp_path, "s1", "qc")
    assert not status.is_completed(tmp_path, "s2", "qc")
    rows = status.get_all(tmp_path)
    assert [r["status"] for r in rows] == ["completed", "failed"]
    assert rows[0]["started_at"] and rows[0]["finished_at"]
