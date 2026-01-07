import pandas as pd
from pathlib import Path


def write_tool_reports(
    parsed_outputs: dict[str, pd.DataFrame], tool_reports: dict[str, Path]
) -> None:
    for tool, df in parsed_outputs.items():
        df.to_csv(tool_reports[tool], sep="\t", index=False)


def write_final_summary(fp: Path, df: pd.DataFrame) -> None:
    df.to_csv(fp, sep="\t", index=False)


def write_assembly_summary(fp: Path, sample_ids: set, sga_version: str) -> None:
    df = pd.DataFrame()
    df["SampleID"] = list(sample_ids)

    try:
        from sunbeam import __version__ as sunbeam_version

        df["sunbeam_version"] = sunbeam_version
    except Exception:
        df["sunbeam_version"] = ""

    if sga_version == "0.0.0":
        df["sbx_sga_version"] = ""
    else:
        df["sbx_sga_version"] = sga_version

    df["metagenomic_sample_id"] = ""
    df["metagenomic_run_id"] = ""
    df["nanopore_path"] = ""
    df["run_number"] = ""
    df["sunbeam_output_path"] = ""
    df["ncbi_id"] = ""

    df.to_csv(fp, sep="\t", index=False)
