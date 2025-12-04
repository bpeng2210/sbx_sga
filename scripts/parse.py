import pandas as pd
import re
from pathlib import Path
from typing import Callable


def _parse_sample_name(fp: Path) -> str:
    return fp.parent.name


def parse_tsv(fp: Path) -> pd.DataFrame:
    df = pd.read_csv(fp, sep="\t")
    df.insert(0, "SampleID", _parse_sample_name(fp))
    return df


def parse_mlst(fp: Path) -> pd.DataFrame:
    # Take in `marc.ast.1076.fa	ecoli_achtman_4	58	adk(6)	fumC(4)	gyrB(4)	icd(16)	mdh(24)	purA(8)	recA(14)`
    # and convert to a dataframe with columns SampleID, classification (e.g. `ecoli_achtman_4 58`), and allele_assignment (e.g. `adk(6) fumC(4) ...`)
    # Note: the tsv comes without a header line
    df = pd.read_csv(fp, sep="\t", header=None)
    df.insert(0, "SampleID", _parse_sample_name(fp))

    classification_row = df.apply(lambda row: f"{row.iloc[2]} {row.iloc[3]}", axis=1)
    df["allele_assignment"] = df.apply(
        lambda row: " ".join(row.iloc[4:].astype(str)), axis=1
    )
    df["classification"] = classification_row

    # Drop unused columns, everything but SampleID, classification, allele_assignment
    df = df[["SampleID", "classification", "allele_assignment"]]

    return df


def parse_bakta_txt(fp: Path) -> pd.DataFrame:
    df = pd.DataFrame()
    with open(fp) as f:
        for line in f:
            if len(line.split(":")) != 2:
                continue
            key, value = line.split(":")
            df[key.strip()] = [value.strip()]

    df.insert(0, "SampleID", _parse_sample_name(fp))
    return df


def _extract_species_name(classification: str) -> str:
    matches = re.findall(r"N[A-Z]_[0-9A-Z]+\.[0-9]", classification)
    if not matches:
        try:
            matches = re.findall(r"[A-Z]{2}_[0-9]+\.[0-9]", classification)
        except Exception as exc:
            print(
                "WARNING: Could not extract species name from classification:",
                classification,
            )
            matches = [" "]

    split_char = matches[0]
    species_split = classification.split(split_char)[1].lstrip()
    species = " ".join(species_split.split()[:2])

    return species


def parse_mash_winning_sorted_tab(
    fp: Path, identity: float, hits: int, median_multiplicity_factor: float
) -> pd.DataFrame:
    df: pd.DataFrame = pd.read_csv(fp, sep="\t", header=None)

    df.columns = [
        "identity",
        "hits_per_thousand",
        "median_multiplicity",
        "val",
        "reference",
        "classification",
    ]

    # Filter by identity threshold
    df = df[df["identity"] >= identity]

    # Filter by hits threshold
    df = df[df["hits_per_thousand"].apply(lambda x: int(x.split("/")[0])) >= hits]

    # Extract species names
    df["species"] = df["classification"].apply(_extract_species_name)

    # Filter out species names that include "phage"
    df = df[~df["species"].str.contains("phage", case=False, na=False)]
    df = df[~df["species"].str.contains("sp.", case=False, na=False)]

    # Filter by median multiplicity factor
    if df.empty:
        return df
    top_median_multiplicity = df["median_multiplicity"].iloc[0]
    df = df[
        df["median_multiplicity"]
        >= top_median_multiplicity * median_multiplicity_factor
    ]

    df.insert(0, "SampleID", _parse_sample_name(fp))
    return df


def _parse_header(h: str) -> dict[str, str]:
    hl = h.split()[1:]  # gets rid of the contig name
    header_dict = dict(item.split("=", 1) for item in hl if "=" in item)
    return header_dict


def parse_fasta(fp: Path) -> pd.DataFrame:
    headers = []
    with open(fp) as f:
        for l in f:
            if l.startswith(">"):
                headers.append(l.strip())

    if not headers:
        return pd.DataFrame(
            {
                "SampleID": [_parse_sample_name(fp)],
                "Total_contigs": [0],
                "Min_coverage": [0],
                "Max_coverage": [0],
                "Total_length": [0],
                "Average_coverage": [0],
            }
        )

    hs = [_parse_header(h) for h in headers]
    total_contigs = len(hs)
    min_cov = min(float(h["cov"]) for h in hs)
    max_cov = max(float(h["cov"]) for h in hs)
    total_length = sum(int(h["len"]) for h in hs)
    avg_cov = sum(int(h["len"]) * float(h["cov"]) for h in hs) / total_length
    rounded_cov = round(avg_cov, 2)

    return pd.DataFrame(
        {
            "SampleID": [_parse_sample_name(fp)],
            "Total_contigs": [total_contigs],
            "Min_coverage": [min_cov],
            "Max_coverage": [max_cov],
            "Total_length": [total_length],
            "Average_coverage": [rounded_cov],
        }
    )


def parse_sylph(fp: Path) -> pd.DataFrame:
    df = parse_tsv(fp)
    df["species"] = df["Contig_name"].apply(_extract_species_name)
    return df


def parse_all_outputs(
    outputs: dict[str, list[Path]],
    parsers: dict[str, Callable],
    parser_kwargs: dict[str, dict] = {},
) -> dict[str, pd.DataFrame]:
    parsed_outputs = {}
    for tool, fps in outputs.items():
        dfs = [parsers[tool](fp, **parser_kwargs.get(tool, {})) for fp in fps]
        combined_df = pd.concat(dfs, ignore_index=True)
        parsed_outputs[tool] = combined_df
    return parsed_outputs
