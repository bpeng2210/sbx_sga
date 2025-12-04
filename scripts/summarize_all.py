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
            from scripts.parse import (
                parse_all_outputs,
                parse_tsv,
                parse_bakta_txt,
                parse_mash_winning_sorted_tab,
                parse_fasta,
                parse_mlst,
                parse_sylph,
            )
            from scripts.write import (
                write_tool_reports,
                write_final_summary,
            )

            parsers = {
                "abritamr": parse_tsv,
                "bakta": parse_bakta_txt,
                "checkm": parse_tsv,
                "mash": parse_mash_winning_sorted_tab,
                "mlst": parse_mlst,
                "shovill": parse_fasta,
                "sylph": parse_sylph,
            }

            outputs: dict[str, list[Path]] = {
                "abritamr": [Path(fp) for fp in snakemake.input.abritamr],  # type: ignore
                "bakta": [Path(fp) for fp in snakemake.input.bakta],  # type: ignore
                "checkm": [Path(fp) for fp in snakemake.input.checkm],  # type: ignore
                "mash": [Path(fp) for fp in snakemake.input.mash],  # type: ignore
                "mlst": [Path(fp) for fp in snakemake.input.mlst],  # type: ignore
                "shovill": [Path(fp) for fp in snakemake.input.shovill],  # type: ignore
                "sylph": [Path(fp) for fp in snakemake.input.sylph],  # type: ignore
            }

            tool_reports = {Path(fp).stem: Path(fp) for fp in snakemake.output.tool_reports}  # type: ignore

            assembly_qcs = snakemake.output.assembly_qcs  # type: ignore
            taxonomic_assignments = snakemake.output.taxonomic_assignments  # type: ignore
            contaminants = snakemake.output.contaminants  # type: ignore
            antimicrobials = snakemake.output.antimicrobials  # type: ignore

            mash_identity = snakemake.params.mash_identity  # type: ignore
            mash_hits = snakemake.params.mash_hits  # type: ignore
            mash_median_multiplicity_factor = snakemake.params.mash_median_multiplicity_factor  # type: ignore

            parser_kwargs = {
                "mash": {
                    "identity": mash_identity,
                    "hits": mash_hits,
                    "median_multiplicity_factor": mash_median_multiplicity_factor,
                }
            }

            # Parse outputs
            log.write("Parsing tool outputs\n")
            parsed_outputs = parse_all_outputs(outputs, parsers, parser_kwargs)

            # Write individual tool reports
            log.write("Writing individual tool reports\n")
            write_tool_reports(parsed_outputs, tool_reports)

            # Produce final summaries
            log.write("Producing final summaries\n")
            write_final_summary(
                assembly_qcs, tools_to_model(parsed_outputs, "assembly_qc")
            )
            write_final_summary(
                taxonomic_assignments,
                tools_to_model(parsed_outputs, "taxonomic_assignment"),
            )
            write_final_summary(
                contaminants, tools_to_model(parsed_outputs, "contaminant")
            )
            write_final_summary(
                antimicrobials, tools_to_model(parsed_outputs, "antimicrobial")
            )
            log.write("Finished writing final summaries\n")
        except Exception as error:
            log.write(f"Encountered error: {error}")
            log.write(traceback.format_exc())
            raise
