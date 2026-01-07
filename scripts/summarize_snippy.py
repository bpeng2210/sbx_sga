import logging
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


if "snakemake" in globals():
    log_fp = snakemake.log[0]  # type: ignore
    with open(log_fp, "w") as log:
        try:
            log.write("Starting summary script\n")
            from scripts.parse import parse_tsv, parse_all_outputs
            from scripts.write import write_tool_reports

            parsers = {
                "snippy": parse_tsv,
            }
            outputs: dict[str, list[Path]] = {
                "snippy": [Path(fp) for fp in snakemake.input.snippy],  # type: ignore
            }

            tool_reports = {"snippy": Path(snakemake.output.snippy)}  # type: ignore

            parsed_outputs = parse_all_outputs(outputs, parsers)
            write_tool_reports(parsed_outputs, tool_reports)
        except Exception:
            log.write("Encountered error during summarization")
            raise
