import pandas as pd
import sys
import traceback
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


if "snakemake" in globals():
    log_fp = snakemake.log[0]  # type: ignore
    with open(log_fp, "w") as log:
        try:
            log.write("Starting summary script\n")
            from scripts.map import tools_to_model
            from scripts.parse import parse_all_outputs, parse_tsv

            parsers = {
                "genomad_plasmid_summary": parse_tsv,
                "genomad_virus_summary": parse_tsv,
                "genomad_plasmid_genes": parse_tsv,
                "genomad_virus_genes": parse_tsv,
            }

            outputs: dict[str, list[Path]] = {
                "genomad_plasmid_summary": [Path(fp) for fp in snakemake.input if "plasmid_summary" in fp],  # type: ignore
                "genomad_virus_summary": [Path(fp) for fp in snakemake.input if "virus_summary" in fp],  # type: ignore
                "genomad_plasmid_genes": [Path(fp) for fp in snakemake.input if "plasmid_genes" in fp],  # type: ignore
                "genomad_virus_genes": [Path(fp) for fp in snakemake.input if "virus_genes" in fp],  # type: ignore
            }

            tool_reports = {Path(fp).stem: Path(fp) for fp in snakemake.output.tool_reports}  # type: ignore

            virus_summary_output = Path(snakemake.output.virus)  # type: ignore

            # Parse outputs
            log.write("Parsing tool outputs\n")
            parsed_outputs = parse_all_outputs(outputs, parsers)

            # Write individual tool reports
            log.write("Writing individual tool reports\n")
            for tool, df in parsed_outputs.items():
                df.to_csv(tool_reports[tool], sep="\t", index=False)

            # Produce final summaries
            log.write(
                "Dummying final virus summary for now (TODO: Implement same as summarize_all.py)\n"
            )
            pd.DataFrame().to_csv(virus_summary_output, sep="\t", index=False)
            log.write("Finished writing final summaries\n")
        except Exception as error:
            log.write(f"Encountered error: {error}")
            log.write(traceback.format_exc())
            raise
