import logging
import pandas as pd
import re
from pathlib import Path
from typing import Callable


logger = logging.getLogger(__name__)


def _parse_sample_name(fp: Path) -> str:
    return fp.parent.name


def parse_tsv(fp: Path) -> pd.DataFrame:
    try:
        logger.debug("Parsing TSV file", extra={"path": str(fp)})
        df = pd.read_csv(fp, sep="\t")
    except pd.errors.EmptyDataError:
        logger.warning("TSV file empty", extra={"path": str(fp)})
        return pd.DataFrame(columns=["SampleID"])

    df.insert(0, "SampleID", _parse_sample_name(fp))
    logger.debug("Parsed TSV dataframe", extra={"path": str(fp)})
    return df


def parse_mlst(fp: Path) -> pd.DataFrame:
    # Take in `marc.ast.1076.fa	ecoli_achtman_4	58	adk(6)	fumC(4)	gyrB(4)	icd(16)	mdh(24)	purA(8)	recA(14)`
    # and convert to a dataframe with columns SampleID, classification (e.g. `ecoli_achtman_4 58`), and allele_assignment (e.g. `adk(6) fumC(4) ...`)
    # Note: the tsv comes without a header line
    try:
        logger.debug("Parsing MLST file", extra={"path": str(fp)})
        df = pd.read_csv(fp, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        logger.warning("MLST file empty", extra={"path": str(fp)})
        return pd.DataFrame(columns=["SampleID", "classification", "allele_assignment"])
    df.insert(0, "SampleID", _parse_sample_name(fp))

    classification_row = df.apply(lambda row: f"{row.iloc[2]} {row.iloc[3]}", axis=1)
    df["allele_assignment"] = df.apply(
        lambda row: " ".join(row.iloc[4:].astype(str)), axis=1
    )
    df["classification"] = classification_row

    # Drop unused columns, everything but SampleID, classification, allele_assignment
    df = df[["SampleID", "classification", "allele_assignment"]]

    logger.debug("Parsed MLST dataframe", extra={"path": str(fp)})
    return df


def parse_bakta_txt(fp: Path) -> pd.DataFrame:
    df = pd.DataFrame()
    logger.debug("Parsing Bakta text file", extra={"path": str(fp)})
    with open(fp) as f:
        for line in f:
            if len(line.split(":")) != 2:
                continue
            key, value = line.split(":")
            df[key.strip()] = [value.strip()]

    df.insert(0, "SampleID", _parse_sample_name(fp))
    logger.debug("Parsed Bakta dataframe", extra={"path": str(fp)})
    return df


def _extract_species_name(classification: str) -> str:
    if not classification:
        return ""

    try:
        cleaned = classification.strip()
        cleaned = re.sub(r"^\[[^\]]*\]\s*", "", cleaned)

        accession_match = re.search(r"[A-Z_]+[0-9]+(?:\.[0-9]+)?", cleaned)
        remainder = (
            cleaned[accession_match.end() :].strip() if accession_match else cleaned
        )

        tokens = remainder.split()
        if len(tokens) >= 2:
            species_tokens = tokens[:2]
        elif tokens:
            species_tokens = [tokens[0]]
        else:
            species_tokens = []

        species = " ".join(token.rstrip(",") for token in species_tokens)
    except Exception as exc:  # pragma: no cover - defensive logging path
        logger.error(
            "Could not extract species name from classification",
            extra={"classification": classification, "exception": str(exc)},
        )
        return ""
    
    return species

def _extract_mash_species(contig_name: str) -> str:
    if not contig_name or pd.isna(contig_name):
        return ""

    extracted = _extract_species_name(contig_name)
    # Filter out classifications we don't want
    if "phage" in species or "sp." in species:
        return ""

def _extract_sylph_species(contig_name: str) -> str:
    if not contig_name or pd.isna(contig_name):
        return ""

    extracted = _extract_species_name(contig_name)
    if not extracted:
        return contig_name
    if "uncultured" in extracted.lower() or "phage" in extracted.lower():
        return contig_name
    return extracted


def parse_mash_winning_sorted_tab(
    fp: Path, identity: float, hits: int, median_multiplicity_factor: float
) -> pd.DataFrame:
    try:
        logger.debug(
            "Parsing mash winning sorted tab",
            extra={
                "path": str(fp),
                "identity_threshold": identity,
                "hits_threshold": hits,
                "median_multiplicity_factor": median_multiplicity_factor,
            },
        )
        df: pd.DataFrame = pd.read_csv(fp, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        logger.warning("Mash file empty", extra={"path": str(fp)})
        return pd.DataFrame(
            columns=[
                "SampleID",
                "identity",
                "hits_per_thousand",
                "median_multiplicity",
                "val",
                "reference",
                "full_classification",
                "species",
            ]
        )

    df.columns = [
        "identity",
        "hits_per_thousand",
        "median_multiplicity",
        "val",
        "reference",
        "full_classification",
    ]

    # Filter by identity threshold
    df = df[df["identity"] >= identity]

    # Filter by hits threshold
    df = df[df["hits_per_thousand"].apply(lambda x: int(x.split("/")[0])) >= hits]

    if df.empty:
        logger.debug("No mash hits after filtering", extra={"path": str(fp)})
        return pd.DataFrame(
            columns=[
                "SampleID",
                "identity",
                "hits_per_thousand",
                "median_multiplicity",
                "val",
                "reference",
                "full_classification",
                "species",
            ]
        )

    # Extract species names
    df["species"] = df["full_classification"].apply(_extract_mash_name)

    # Filter by median multiplicity factor
    if df.empty:
        logger.debug("No mash hits after species filtering", extra={"path": str(fp)})
        return pd.DataFrame(
            columns=[
                "SampleID",
                "identity",
                "hits_per_thousand",
                "median_multiplicity",
                "val",
                "reference",
                "full_classification",
                "species",
            ]
        )
    top_median_multiplicity = df["median_multiplicity"].iloc[0]
    df = df[
        df["median_multiplicity"]
        >= top_median_multiplicity * median_multiplicity_factor
    ]

    # Sort by hits with the largest at the top
    df = df.sort_values(
        by="hits_per_thousand",
        key=lambda col: col.apply(lambda x: int(x.split("/")[0])),
        ascending=False,
    )
    # Remove duplicate species assignments, keeping the top hit
    df = df.drop_duplicates(subset=["species"], keep="first")

    df.insert(0, "SampleID", _parse_sample_name(fp))
    logger.debug("Parsed mash dataframe", extra={"path": str(fp)})
    return df


def _parse_header(h: str) -> dict[str, str]:
    hl = h.split()[1:]  # gets rid of the contig name
    header_dict = dict(item.split("=", 1) for item in hl if "=" in item)
    return header_dict


def parse_fasta(fp: Path) -> pd.DataFrame:
    logger.debug("Parsing FASTA file", extra={"path": str(fp)})
    headers = []
    with open(fp) as f:
        for l in f:
            if l.startswith(">"):
                headers.append(l.strip())

    if not headers:
        logger.debug("FASTA has no headers", extra={"path": str(fp)})
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

    df = pd.DataFrame(
        {
            "SampleID": [_parse_sample_name(fp)],
            "Total_contigs": [total_contigs],
            "Min_coverage": [min_cov],
            "Max_coverage": [max_cov],
            "Total_length": [total_length],
            "Average_coverage": [rounded_cov],
        }
    )

    logger.debug("Parsed FASTA dataframe", extra={"path": str(fp)})
    return df


def parse_sylph(fp: Path) -> pd.DataFrame:
    df = parse_tsv(fp)

    if "Contig_name" not in df.columns:
        df["Contig_name"] = pd.Series(dtype=str)

    if df.empty:
        df["species"] = pd.Series(dtype=str)
        return df[["SampleID", "Contig_name", "species"]]

    df["species"] = df["Contig_name"].apply(_extract_sylph_species)
    logger.debug("Parsed sylph dataframe", extra={"path": str(fp)})
    return df


def parse_all_outputs(
    outputs: dict[str, list[Path]],
    parsers: dict[str, Callable],
    parser_kwargs: dict[str, dict] = {},
) -> dict[str, pd.DataFrame]:
    parsed_outputs = {}
    for tool, fps in outputs.items():
        logger.debug(
            "Parsing outputs for tool",
            extra={"tool": tool, "file_paths": [str(fp) for fp in fps]},
        )
        dfs = []
        for fp in fps:
            try:
                dfs.append(parsers[tool](fp, **parser_kwargs.get(tool, {})))
            except Exception:
                logger.exception(
                    "Failed to parse file",
                    extra={"tool": tool, "path": str(fp)},
                )
                continue

        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
        else:
            logger.error(
                "No successful parses for tool",
                extra={"tool": tool, "file_paths": [str(fp) for fp in fps]},
            )
            combined_df = pd.DataFrame()

        logger.debug("Combined dataframe for tool", extra={"tool": tool})
        parsed_outputs[tool] = combined_df
    return parsed_outputs
