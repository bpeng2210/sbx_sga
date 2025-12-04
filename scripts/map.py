import pandas as pd


def _merge_dfs_on_sample_id(dfs: list[pd.DataFrame]) -> pd.DataFrame:
    """Merge dataframes on SampleID while coalescing overlapping columns.

    When multiple tools map to the same target fields (e.g., ``classification``),
    pandas' default merge would create suffixed columns (``_x``/``_y``). This
    helper performs an outer merge and then combines overlapping columns so the
    resulting dataframe retains a single shared column for each field.
    """

    dfs = list(dfs)
    if not dfs:
        return pd.DataFrame(columns=["SampleID"])

    merged = dfs[0]
    for df in dfs[1:]:
        overlapping_cols = [
            col for col in merged.columns if col in df.columns and col != "SampleID"
        ]
        merged = pd.merge(merged, df, on="SampleID", how="outer", suffixes=("", "_dup"))

        for col in overlapping_cols:
            dup_col = f"{col}_dup"
            merged[col] = merged[col].combine_first(merged[dup_col])
            merged = merged.drop(columns=[dup_col])

        dup_cols = [col for col in merged.columns if col.endswith("_dup")]
        if dup_cols:
            merged = merged.drop(columns=dup_cols)

    return merged


def tools_to_assembly_qc(parsed_outputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    shovill_df = parsed_outputs.get("shovill", pd.DataFrame())
    shovill_df = shovill_df.rename(
        columns={
            "Total_contigs": "contig_count",
            "Min_coverage": "min_contig_coverage",
            "Max_coverage": "max_contig_coverage",
            "Total_length": "genome_size",
            "Average_coverage": "avg_contig_coverage",
        }
    )

    bakta_df = parsed_outputs.get("bakta", pd.DataFrame())
    bakta_df = bakta_df.rename(
        columns={"GC": "gc_content", "N50": "n50", "CDSs": "cds"}
    )

    checkm_df = parsed_outputs.get("checkm", pd.DataFrame())
    checkm_df = checkm_df.rename(
        columns={"Completeness": "completeness", "Contamination": "contamination"}
    )

    df = _merge_dfs_on_sample_id(
        [df for df in [shovill_df, bakta_df, checkm_df] if not df.empty]
    )
    return df[
        [
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
    ]


def tools_to_taxonomic_assignment(
    parsed_outputs: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    mlst_df = parsed_outputs.get("mlst", pd.DataFrame())
    mlst_df = mlst_df.rename(
        columns={
            "classification": "classification",
            "allele_assignment": "comment",
        }
    )
    # Add tool field
    if not mlst_df.empty:
        mlst_df["tool"] = "mlst"

    sylph_df = parsed_outputs.get("sylph", pd.DataFrame())
    sylph_df = sylph_df.rename(
        columns={
            "Contig_name": "classification",
        }
    )
    # Add empty comment column for sylph if it doesn't exist
    if not sylph_df.empty:
        if "comment" not in sylph_df.columns:
            sylph_df["comment"] = ""
        sylph_df["tool"] = "sylph"

    # Concatenate instead of merge to allow multiple classifications per sample
    dfs_to_concat = [df for df in [mlst_df, sylph_df] if not df.empty]
    if not dfs_to_concat:
        return pd.DataFrame(columns=["SampleID"])

    return pd.concat(dfs_to_concat, ignore_index=True)[
        ["SampleID", "classification", "comment", "tool"]
    ]


def tools_to_contaminant(parsed_outputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    mash_df = parsed_outputs.get("mash", pd.DataFrame())
    mash_df = mash_df.rename(
        columns={
            "hits_per_thousand": "confidence",
            "species": "classification",
        }
    )
    # Add tool field
    if not mash_df.empty:
        mash_df["tool"] = "mash"

    # Use concatenation approach like taxonomic assignment for consistency
    dfs_to_concat = [df for df in [mash_df] if not df.empty]
    if not dfs_to_concat:
        return pd.DataFrame(columns=["SampleID"])

    return pd.concat(dfs_to_concat, ignore_index=True)[
        ["SampleID", "classification", "confidence", "tool"]
    ]


def tools_to_antimicrobial(parsed_outputs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    abritamr_df = parsed_outputs.get("abritamr", pd.DataFrame())
    abritamr_df = abritamr_df.rename(
        columns={
            "Contig id": "contig_id",
            "Gene symbol": "gene_symbol",
            "Sequence name": "gene_name",
            "Accession of closest sequence": "accession",
            "Element type": "element_type",
            "Subclass": "resistance_product",
        }
    )

    df = _merge_dfs_on_sample_id([df for df in [abritamr_df] if not df.empty])
    return df[
        [
            "SampleID",
            "contig_id",
            "gene_symbol",
            "gene_name",
            "accession",
            "element_type",
            "resistance_product",
        ]
    ]


def tools_to_model(parsed_outputs: dict[str, pd.DataFrame], model: str) -> pd.DataFrame:
    if model == "assembly_qc":
        return tools_to_assembly_qc(parsed_outputs)
    elif model == "taxonomic_assignment":
        return tools_to_taxonomic_assignment(parsed_outputs)
    elif model == "contaminant":
        return tools_to_contaminant(parsed_outputs)
    elif model == "antimicrobial":
        return tools_to_antimicrobial(parsed_outputs)
    else:
        raise ValueError(f"Unknown model: {model}")
