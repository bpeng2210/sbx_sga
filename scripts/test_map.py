import pandas as pd
from scripts.map import tools_to_model


parsed_outputs = {
    "shovill": pd.DataFrame(
        [
            {
                "SampleID": "S1",
                "Total_contigs": 10,
                "Min_coverage": 5,
                "Max_coverage": 40,
                "Total_length": 50000,
                "Average_coverage": 25,
            },
            {
                "SampleID": "S2",
                "Total_contigs": 20,
                "Min_coverage": 4,
                "Max_coverage": 35,
                "Total_length": 75000,
                "Average_coverage": 30,
            },
        ]
    ),
    "bakta": pd.DataFrame(
        [
            {"SampleID": "S1", "GC": 51.1, "N50": 1500, "CDSs": "2000"},
            {"SampleID": "S3", "GC": 47.0, "N50": 900, "CDSs": "1800"},
        ]
    ),
    "checkm": pd.DataFrame(
        [
            {"SampleID": "S2", "Completeness": 98.5, "Contamination": 1.5},
        ]
    ),
    "mash": pd.DataFrame(
        [
            {
                "SampleID": "S1",
                "hits_per_thousand": "911/1000",
                "species": "",
            },
            {
                "SampleID": "S3",
                "hits_per_thousand": "110/1000",
                "species": "Bug",
            },
        ]
    ),
    "mlst": pd.DataFrame(
        [
            {
                "SampleID": "S1",
                "classification": "schema1_1 10",
                "allele_assignment": "geneA(1) geneB(2)",
            },
            {
                "SampleID": "S3",
                "classification": "schema2_2 20",
                "allele_assignment": "geneC(3) geneD(4)",
            },
        ]
    ),
    "sylph": pd.DataFrame(
        [
            {
                "SampleID": "S2",
                "species": "schema_st_15",
            },
            {
                "SampleID": "S3",
                "species": "schema_st_25",
            },
        ]
    ),
    "abritamr": pd.DataFrame(
        [
            {
                "SampleID": "S2",
                "Contig id": "contig_1",
                "Gene symbol": "geneA",
                "Sequence name": "geneA_full",
                "Accession of closest sequence": "acc1",
                "Element type": "element",
                "Subclass": "resistance_product",
            }
        ]
    ),
}


def test_tools_to_assembly_qc():
    assembly_qc_df = tools_to_model(parsed_outputs, "assembly_qc")

    assert not assembly_qc_df.empty
    assert set(assembly_qc_df["SampleID"].tolist()) == {"S1", "S2", "S3"}
    assert all(
        c in assembly_qc_df.columns for c in ["contig_count", "gc_content", "n50"]
    )


def test_tools_to_taxonomic_assignment():
    taxonomic_assignment_df = tools_to_model(parsed_outputs, "taxonomic_assignment")

    assert not taxonomic_assignment_df.empty
    assert set(taxonomic_assignment_df["SampleID"].tolist()) == {"S1", "S2", "S3"}
    print("HERE:", taxonomic_assignment_df)
    assert all(
        c in taxonomic_assignment_df.columns
        for c in ["classification", "comment", "tool"]
    )
    # Get both entries for S3
    s3_entries = taxonomic_assignment_df[taxonomic_assignment_df["SampleID"] == "S3"]
    assert len(s3_entries) == 2


def test_tools_to_contaminant():
    contaminant_df = tools_to_model(parsed_outputs, "contaminant")

    assert not contaminant_df.empty
    assert set(contaminant_df["SampleID"].tolist()) == {"S1", "S3"}
    assert all(
        c in contaminant_df.columns for c in ["classification", "confidence", "tool"]
    )


def test_tools_to_antimicrobial():
    antimicrobial_df = tools_to_model(parsed_outputs, "antimicrobial")

    assert not antimicrobial_df.empty
    assert set(antimicrobial_df["SampleID"].tolist()) == {"S2"}
    print(antimicrobial_df)
    assert all(
        [
            c in antimicrobial_df.columns
            for c in [
                "contig_id",
                "gene_symbol",
                "gene_name",
                "accession",
                "element_type",
                "resistance_product",
            ]
        ]
    )


def test_tools_handle_empty_inputs():
    empty_parsed_outputs = {
        "shovill": pd.DataFrame(
            columns=[
                "SampleID",
                "Total_contigs",
                "Min_coverage",
                "Max_coverage",
                "Total_length",
                "Average_coverage",
            ]
        ),
        "bakta": pd.DataFrame(columns=["SampleID", "GC", "N50", "CDSs"]),
        "checkm": pd.DataFrame(columns=["SampleID", "Completeness", "Contamination"]),
        "mash": pd.DataFrame(columns=["SampleID", "species", "hits_per_thousand"]),
        "mlst": pd.DataFrame(
            columns=["SampleID", "classification", "allele_assignment"]
        ),
        "sylph": pd.DataFrame(columns=["SampleID", "species"]),
        "abritamr": pd.DataFrame(
            columns=[
                "SampleID",
                "Contig id",
                "Gene symbol",
                "Sequence name",
                "Accession of closest sequence",
                "Element type",
                "Subclass",
            ]
        ),
    }

    assembly_qc_df = tools_to_model(empty_parsed_outputs, "assembly_qc")
    assert assembly_qc_df.empty
    assert list(assembly_qc_df.columns) == [
        "SampleID",
        "contig_count",
        "min_contig_coverage",
        "max_contig_coverage",
        "genome_size",
        "avg_contig_coverage",
        "gc_content",
        "n50",
        "cds",
        "completeness",
        "contamination",
    ]

    taxonomic_assignment_df = tools_to_model(
        empty_parsed_outputs, "taxonomic_assignment"
    )
    assert taxonomic_assignment_df.empty
    assert list(taxonomic_assignment_df.columns) == [
        "SampleID",
        "classification",
        "comment",
        "tool",
    ]

    contaminant_df = tools_to_model(empty_parsed_outputs, "contaminant")
    assert contaminant_df.empty
    assert list(contaminant_df.columns) == [
        "SampleID",
        "classification",
        "confidence",
        "tool",
    ]

    antimicrobial_df = tools_to_model(empty_parsed_outputs, "antimicrobial")
    assert antimicrobial_df.empty
    assert list(antimicrobial_df.columns) == [
        "SampleID",
        "contig_id",
        "gene_symbol",
        "gene_name",
        "accession",
        "element_type",
        "resistance_product",
    ]
