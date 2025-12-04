import pandas as pd
from pathlib import Path


def write_tool_reports(
    parsed_outputs: dict[str, pd.DataFrame], tool_reports: dict[str, Path]
) -> None:
    for tool, df in parsed_outputs.items():
        df.to_csv(tool_reports[tool], sep="\t", index=False)


def write_final_summary(fp: Path, df: pd.DataFrame) -> None:
    df.to_csv(fp, sep="\t", index=False)
