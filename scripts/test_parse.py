import pytest
import shutil
from pathlib import Path
from .parse import (
    parse_tsv,
    parse_bakta_txt,
    parse_mash_winning_sorted_tab,
    parse_fasta,
    parse_mlst,
    parse_sylph,
    _extract_species_name,
)


@pytest.fixture
def test_reports_fp():
    return Path(__file__).parent.parent / ".tests/data/test_reports"


@pytest.fixture
def sample_report_fp(tmp_path, test_reports_fp):
    def _factory(tool: str, sample: str, filename: str) -> Path:
        dest_dir = tmp_path / tool / sample
        dest_dir.mkdir(parents=True, exist_ok=True)
        src = test_reports_fp / tool / filename
        dest = dest_dir / filename
        shutil.copy(src, dest)
        return dest

    return _factory


### AbritAMR ###
def test_amr(sample_report_fp):
    sample_name = "dummy"
    fp = sample_report_fp("abritamr", sample_name, "amrfinder.out")
    df = parse_tsv(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "Contig id" in df.columns
    assert "Gene symbol" in df.columns
    assert "Subclass" in df.columns


def test_parse_tsv_empty(tmp_path):
    sample_name = "empty"
    fp = tmp_path / "abritamr" / sample_name / "amrfinder.out"
    fp.parent.mkdir(parents=True, exist_ok=True)
    fp.write_text("")

    df = parse_tsv(fp)

    assert df.empty
    assert list(df.columns) == ["SampleID"]


def test_parse_tsv_header_only(tmp_path):
    sample_name = "header_only"
    fp = tmp_path / "abritamr" / sample_name / "amrfinder.out"
    fp.parent.mkdir(parents=True, exist_ok=True)
    fp.write_text("ColA\tColB\n")

    df = parse_tsv(fp)

    assert df.empty
    assert list(df.columns) == ["SampleID", "ColA", "ColB"]


### Bakta ###
def test_bakta(sample_report_fp):
    sample_name = "marc.ast.1076"
    fp = sample_report_fp("bakta", sample_name, "marc.ast.1076.txt")
    df = parse_bakta_txt(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "CDSs" in df.columns
    assert df["SampleID"].unique().tolist() == [sample_name]
    assert df["CDSs"].iloc[0] == "4462"


### CheckM ###
def test_checkm(sample_report_fp):
    sample_name = "dummy"
    fp = sample_report_fp("checkm", sample_name, "quality_report.tsv")
    df = parse_tsv(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "Completeness" in df.columns
    assert "Contamination" in df.columns


### Mash ###
def test_parse_mash_marc_3111(sample_report_fp):
    sample_name = "marc.bacteremia.3111"
    fp = sample_report_fp(
        "mash",
        sample_name,
        "marc.bacteremia.3111_sorted_winning.tab",
    )
    df = parse_mash_winning_sorted_tab(
        fp, identity=0.85, hits=100, median_multiplicity_factor=0.05
    )
    assert not df.empty
    assert "species" in df.columns
    assert all(df["identity"] >= 0.85)
    assert all(df["hits_per_thousand"].apply(lambda x: int(x.split("/")[0])) >= 100)
    assert df["SampleID"].unique().tolist() == [sample_name]
    assert df["species"].iloc[0] == "Enterobacter cloacae"
    assert df["species"].iloc[1] == "Enterobacter kobei"


def test_parse_mash_s234_ori(sample_report_fp):
    sample_name = "s234.ori.lightblue.b"
    fp = sample_report_fp(
        "mash",
        sample_name,
        "s234.ori.lightblue.b_sorted_winning.tab",
    )
    df = parse_mash_winning_sorted_tab(
        fp, identity=0.85, hits=100, median_multiplicity_factor=0.05
    )
    assert not df.empty
    assert all(df["identity"] >= 0.85)
    assert all(df["hits_per_thousand"].apply(lambda x: int(x.split("/")[0])) >= 100)
    assert df["SampleID"].unique().tolist() == [sample_name]
    assert df["species"].iloc[0] == "Bacillus cereus"
    assert df["species"].iloc[1] == "Pseudomonas denitrificans"
    assert df["species"].iloc[2] == "Stenotrophomonas maltophilia"
    assert df["species"].iloc[3] == "Stenotrophomonas acidaminiphila"


def test_parse_mash_marc_235(sample_report_fp):
    sample_name = "marc.entero.235"
    fp = sample_report_fp(
        "mash",
        sample_name,
        "marc.entero.235_sorted_winning.tab",
    )
    df = parse_mash_winning_sorted_tab(
        fp, identity=0.85, hits=100, median_multiplicity_factor=0.05
    )
    assert df.empty

    df = parse_mash_winning_sorted_tab(
        fp, identity=0.80, hits=10, median_multiplicity_factor=0.05
    )
    assert not df.empty
    assert all(df["identity"] >= 0.80)
    assert all(df["hits_per_thousand"].apply(lambda x: int(x.split("/")[0])) >= 10)
    assert df["SampleID"].unique().tolist() == [sample_name]
    assert df["species"].iloc[0] == "Serratia marcescens"


def test_parse_mash_empty(tmp_path):
    sample_name = "empty"
    fp = tmp_path / "mash" / sample_name / "sample.tab"
    fp.parent.mkdir(parents=True, exist_ok=True)
    fp.write_text("")

    df = parse_mash_winning_sorted_tab(
        fp, identity=0.85, hits=100, median_multiplicity_factor=0.05
    )

    assert df.empty
    assert "SampleID" in df.columns


### MLST ###
def test_mlst_1076(sample_report_fp):
    sample_name = "marc.ast.1076"
    fp = sample_report_fp("mlst", sample_name, "marc.ast.1076.mlst")
    df = parse_mlst(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "classification" in df.columns
    assert "allele_assignment" in df.columns
    assert df["SampleID"].unique().tolist() == [sample_name]
    assert df["classification"].iloc[0] == "ecoli_achtman_4 58"
    assert (
        df["allele_assignment"].iloc[0]
        == "adk(6) fumC(4) gyrB(4) icd(16) mdh(24) purA(8) recA(14)"
    )


def test_mlst_empty(tmp_path):
    sample_name = "empty"
    fp = tmp_path / "mlst" / sample_name / "sample.mlst"
    fp.parent.mkdir(parents=True, exist_ok=True)
    fp.write_text("")

    df = parse_mlst(fp)

    assert df.empty
    assert list(df.columns) == ["SampleID", "classification", "allele_assignment"]


### Shovill ###
def test_fasta(sample_report_fp):
    sample_name = "marc.ast.1076"
    fp = sample_report_fp("shovill", sample_name, "marc.ast.1076.fa")
    df = parse_fasta(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "Total_contigs" in df.columns
    assert df["SampleID"].unique().tolist() == [sample_name]


def test_empty_fasta(sample_report_fp):
    sample_name = "empty"
    fp = sample_report_fp("shovill", sample_name, "empty.fa")
    df = parse_fasta(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert df["SampleID"].unique().tolist() == [sample_name]
    assert df["Total_contigs"].iloc[0] == 0


### Sylph ###
def test_parse_tsv_sylph_3151(sample_report_fp):
    sample_name = "marc.bacteremia.3151"
    fp = sample_report_fp("sylph", sample_name, "marc.bacteremia.3151.tsv")
    df = parse_tsv(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "Contig_name" in df.columns
    assert df["SampleID"].unique().tolist() == [sample_name]
    assert "Actinomyces" in df["Contig_name"].iloc[0]


def test_parse_tsv_sylph_s234_ori(sample_report_fp):
    sample_name = "s234.ori.lightblue.b"
    fp = sample_report_fp("sylph", sample_name, "s234.ori.lightblue.b.tsv")
    df = parse_tsv(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "Contig_name" in df.columns
    assert df["SampleID"].unique().tolist() == [sample_name]
    assert df.shape[0] > 5


def test_parse_sylph_empty(tmp_path):
    sample_name = "empty"
    fp = tmp_path / "sylph" / sample_name / "sample.tsv"
    fp.parent.mkdir(parents=True, exist_ok=True)
    fp.write_text("")

    df = parse_sylph(fp)

    assert df.empty
    assert list(df.columns) == ["SampleID", "Contig_name", "species"]


def test_extract_species_name_variations():
    cases = [
        (
            "BMZV01000001.1 Kocuria marina KCTC 9943 DNA, sequence01, whole genome shotgun sequence",
            "Kocuria marina",
        ),
        (
            "CP007494.1 Bordetella holmesii ATCC 51541, complete genome",
            "Bordetella holmesii",
        ),
        (
            "NZ_SCMK01000010.1 Sphingomonas sp. 1F27F7B NODE_10_length_165275_cov_23.714959, whole genome shotgun sequence",
            "Sphingomonas sp.",
        ),
        (
            "[4 seqs] NC_013922.1 Natrialba magadii ATCC 43099, complete genome [...]",
            "Natrialba magadii",
        ),
        (
            "NC_013947.1 Stackebrandtia nassauensis DSM 44728, complete genome",
            "Stackebrandtia nassauensis",
        ),
    ]

    for classification, expected in cases:
        assert _extract_species_name(classification) == expected


### Snippy ###
def test_snippy(sample_report_fp):
    sample_name = "sample1"
    fp = sample_report_fp("snippy", sample_name, "snps.tab")
    df = parse_tsv(fp)

    assert not df.empty
    assert "SampleID" in df.columns
    assert "CHROM" in df.columns
    assert "POS" in df.columns


def test_snippy_empty(sample_report_fp):
    sample_name = "sample2"
    fp = sample_report_fp("snippy", sample_name, "empty.tab")
    df = parse_tsv(fp)

    assert df.empty
