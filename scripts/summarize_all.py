import pandas as pd
import traceback
from pathlib import Path


if "snakemake" in globals():
    log_fp = snakemake.log[0]  # type: ignore
    with open(log_fp, "w") as log:
        try:
            log.write("Starting summary script\n")
            from .map import antimicrobial, assembly_qc, reduce_dataframe, taxonomic_assignment
            from .parse import (
                parse_all_outputs,
                parse_tsv,
                parse_bakta_txt,
                parse_mash_winning_sorted_tab,
                parse_fasta,
            )
            
            parsers = {
                "abritamr": parse_tsv,
                "bakta": parse_bakta_txt,
                "checkm": parse_tsv,
                "mash": parse_mash_winning_sorted_tab,
                "mlst": parse_tsv,
                "shovill": parse_fasta,
                "sylph": parse_tsv,
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
            if set(tool_reports.keys()) != set(parsed_outputs.keys()):
                log.write(
                    f"Warning: tool reports keys {list(tool_reports.keys())} do not match parsed outputs keys {list(parsed_outputs.keys())}\n"
                )
            for tool, df in parsed_outputs.items():
                df.to_csv(tool_reports[tool], sep="\t", index=False)

            # Produce final summaries
            log.write("Producing final summaries\n")
            assembly_qc_df = pd.merge(
                *[
                    reduce_dataframe(df, tool)
                    for tool, df in parsed_outputs.items()
                    if tool in assembly_qc
                ],
                on="SampleID",
                how="outer",
            )
            assembly_qc_df.to_csv(assembly_qcs, sep="\t", index=False)

            taxonomic_assignment_df = pd.merge(
                *[
                    reduce_dataframe(df, tool)
                    for tool, df in parsed_outputs.items()
                    if tool in taxonomic_assignment
                ],
                on="SampleID",
                how="outer",
            )
            taxonomic_assignment_df.to_csv(taxonomic_assignments, sep="\t", index=False)

            antimicrobial_df = pd.merge(
                *[
                    reduce_dataframe(df, tool)
                    for tool, df in parsed_outputs.items()
                    if tool in antimicrobial
                ],
                on="SampleID",
                how="outer",
            )
            antimicrobial_df.to_csv(antimicrobials, sep="\t", index=False)
            log.write("Finished writing final summaries\n")
        except Exception as error:
            log.write(f"Encountered error: {error}")
            log.write(traceback.format_exc())
            raise
