"""Tests for CLI helpers: sample-ID inference, read-pair discovery, commands."""
from pathlib import Path

import pytest
from click.testing import CliRunner

from phageflow import __version__
from phageflow.cli import (
    _find_read_pairs,
    _infer_sample_id,
    _infer_sample_id_from_fasta,
    _r2_name_for,
    _sample_id_from_genome_path,
    cli,
)


@pytest.mark.parametrize("r1, expected", [
    ("sampleA_R1_001.fastq.gz", "sampleA"),
    ("sampleA_R1.fastq.gz", "sampleA"),
    ("sampleA_1.fastq.gz", "sampleA"),
    ("sampleA_1.fq.gz", "sampleA"),
    ("SRR22410879_R1.fastq.gz", "SRR22410879"),
    ("weird.fastq", "weird"),
])
def test_infer_sample_id(r1, expected):
    assert _infer_sample_id(r1) == expected


@pytest.mark.parametrize("fasta, expected", [
    ("SRR22410879_contigs.fasta", "SRR22410879"),
    ("SRR22410879_contigs.fasta_virus.fna", "SRR22410879"),
    ("SRR22410879_virus.fna", "SRR22410879"),
    ("sample_contigs.fa", "sample"),
])
def test_infer_sample_id_from_fasta(fasta, expected):
    assert _infer_sample_id_from_fasta(fasta) == expected


@pytest.mark.parametrize("r1, r2", [
    # The "R1" inside an SRA accession must not be rewritten.
    ("SRR1234_R1.fastq.gz", "SRR1234_R2.fastq.gz"),
    ("ERR1001_R1.fq.gz", "ERR1001_R2.fq.gz"),
    ("S1_L001_R1_001.fastq.gz", "S1_L001_R2_001.fastq.gz"),
    ("sample_1.fastq.gz", "sample_2.fastq.gz"),
    ("lib10_R1.fastq", "lib10_R2.fastq"),
    ("lane_1_R1.fastq", "lane_1_R2.fastq"),
    ("x.read1.fastq", "x.read2.fastq"),
])
def test_r2_name_for(r1, r2):
    assert _r2_name_for(r1) == r2


@pytest.mark.parametrize("name", [
    "sample_R2.fastq.gz", "sample_2.fastq.gz", "lane_1_R2.fastq",
    "sample.fastq.gz", "notes.txt",
])
def test_r2_name_for_non_r1(name):
    assert _r2_name_for(name) is None


def test_find_read_pairs(tmp_path: Path):
    for name in ("SRR1234_R1.fastq.gz", "SRR1234_R2.fastq.gz",
                 "B_1.fq.gz", "B_2.fq.gz", "orphan_R1.fastq.gz", "README.md"):
        (tmp_path / name).write_text("x")
    pairs = _find_read_pairs(tmp_path)
    names = sorted((a.name, b.name) for a, b in pairs)
    assert names == [
        ("B_1.fq.gz", "B_2.fq.gz"),
        ("SRR1234_R1.fastq.gz", "SRR1234_R2.fastq.gz"),
    ]


def test_sample_id_from_genome_path(tmp_path: Path):
    g = tmp_path / "SRR1" / "annotation_ready" / "phages" / "Phage_candidate_001.fasta"
    assert _sample_id_from_genome_path(str(g)) == "SRR1"
    assert _sample_id_from_genome_path(str(tmp_path / "x.fasta")) == "x"


def test_version_and_help():
    runner = CliRunner()
    res = runner.invoke(cli, ["--version"])
    assert res.exit_code == 0 and __version__ in res.output
    res = runner.invoke(cli, ["--help"])
    assert res.exit_code == 0
    # Commands whose modules do not exist must not be exposed.
    assert "safety" not in res.output.split("Commands:")[-1]


def test_run_from_module_choices_match_runner():
    from phageflow.modules.runner import _MODULE_ORDER
    run_cmd = cli.commands["run"]
    opt = next(p for p in run_cmd.params if p.name == "from_module")
    assert list(opt.type.choices) == _MODULE_ORDER


def test_output_dir_required_without_project(tmp_path: Path):
    r1 = tmp_path / "s_R1.fastq"; r1.write_text("@a\nA\n+\nI\n")
    r2 = tmp_path / "s_R2.fastq"; r2.write_text("@a\nA\n+\nI\n")
    res = CliRunner().invoke(cli, ["qc", "--r1", str(r1), "--r2", str(r2)])
    assert res.exit_code == 1


def test_init_and_status(tmp_path: Path):
    runner = CliRunner()
    proj = tmp_path / "proj"
    res = runner.invoke(cli, ["init", str(proj)])
    assert res.exit_code == 0
    assert (proj / "config" / "config.yaml").exists()
    for sub in ("raw", "results", "reports"):
        assert (proj / sub).is_dir()
    res = runner.invoke(cli, ["status", "--project", str(proj)])
    assert res.exit_code == 0
