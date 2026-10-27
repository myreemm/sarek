#!/usr/bin/env python3

import argparse
import os
import sys
import traceback

## AlphaGenome wrapper imports

from alphagenome_wrapper.result_io import (
    write_summary_scores,
    write_tidy_scores,
)

from alphagenome_wrapper.vcf_io import (
    InfoField,
    open_annotated_vcf,
    set_info_value,
)

from alphagenome_wrapper.scoring import (
    AlphaGenomeScoreResult,
    compute_alphagenome_score,
    get_selected_scorers,
)

from alphagenome_wrapper.variants import (
    VariantClass,
    VariantRecordData,
    extract_variant_record_data,
    is_scorable_variant_class,
)
from alphagenome_wrapper.aggregation import ( 
    CHROMATIN_GROUP,
    EXPRESSION_GROUP,
    empty_score_summary, )

from dataclasses import dataclass
from enum import Enum

# AlphaGenome specific imports
from alphagenome.models import dna_client

# Constants:

# Missing value used in VCF and TSV outputs
MISSING_VALUE = "."

# Enums definitions:

class VariantStatus(str, Enum):
    """Processing states assigned by the AlphaGenome wrapper."""

    SCORED = "SCORED"
    ERROR = "ERROR"
    SKIPPED_UNSUPPORTED = "SKIPPED_UNSUPPORTED_VARIANT"
    SKIPPED_LIMIT = "SKIPPED_MAX_VARIANTS"


# Dataclass:


@dataclass(frozen=True)
class VariantOutputValues:
    """Values written to the annotated VCF and summary TSV."""

    # Global top differential result
    max_abs_score: float | str
    top_score: float | str
    top_output_type: str
    top_track: str
    top_gene: str

    # Chromatin-accessibility result
    chromatin_score: float | str
    chromatin_abs_score: float | str
    chromatin_output_type: str
    chromatin_track: str

    # Gene-expression result
    expression_score: float | str
    expression_abs_score: float | str
    expression_gene: str
    expression_track: str


@dataclass(frozen=True)
class ProcessedVariantResult:
    """Result of processing one VCF variant record."""

    record: VariantRecordData
    score_result: AlphaGenomeScoreResult | None
    status: VariantStatus


 
# Load AlphaGenome API key from file or environment variable.
# If an API key file is explicitly provided, it must exist and must not be empty
def load_api_key(api_key_file):
    if api_key_file:
        if not os.path.exists(api_key_file):
            raise FileNotFoundError(
                f"AlphaGenome API key file does not exist: {api_key_file}"
            )

        with open(api_key_file, "r") as handle:
            api_key = handle.read().strip()

        if not api_key:
            raise RuntimeError(
                f"AlphaGenome API key file is empty: {api_key_file}"
            )

        if api_key.lower() in {"dummy_key", "your_api_key_here", "changeme"}:
            raise RuntimeError(
                f"AlphaGenome API key file contains a placeholder value: {api_key_file}"
            )

        return api_key

    api_key = os.environ.get("ALPHA_GENOME_API_KEY", "").strip()

    if not api_key:
        raise RuntimeError(
            "No AlphaGenome API key provided. "
            "Use --api-key-file or set ALPHA_GENOME_API_KEY."
        )

    return api_key


# parse configuration (organism & sequence length) and resolve to AlphaGenome API values
def resolve_organism(organism_str):
    organism_map = {
        "human": dna_client.Organism.HOMO_SAPIENS,
        "mouse": dna_client.Organism.MUS_MUSCULUS,
    }

    key = organism_str.lower()

    if key not in organism_map:
        raise ValueError(f"Unsupported organism: {organism_str}")

    return organism_map[key]

def resolve_sequence_length(sequence_length_str):
    seq_len_key = f"SEQUENCE_LENGTH_{sequence_length_str}"

    if seq_len_key not in dna_client.SUPPORTED_SEQUENCE_LENGTHS:
        raise ValueError(f"Unsupported sequence length: {sequence_length_str}")

    return dna_client.SUPPORTED_SEQUENCE_LENGTHS[seq_len_key]


def process_variant_record(
    variant,
    *,
    dna_model,
    organism,
    sequence_length,
    selected_scorers,
    scored_count: int,
    max_variants: int | None,
) -> ProcessedVariantResult:
    """
    Process one VCF record and optionally perform AlphaGenome scoring.

    The function classifies the variant, checks whether it is supported,
    enforces the configured variant limit, and performs AlphaGenome scoring
    when applicable.

    Args:
        variant:
            Variant record parsed by cyvcf2.
        dna_model:
            Initialized AlphaGenome API client.
        organism:
            AlphaGenome organism configuration.
        sequence_length:
            Sequence context length used for prediction.
        selected_scorers:
            AlphaGenome variant scorers selected for the current run.
        scored_count:
            Number of successfully scored variants before this record.
        max_variants:
            Maximum number of variants that may be scored. ``None`` means
            that no explicit limit is applied.

    Returns:
        Structured processing result for the current variant.
    """

    record = extract_variant_record_data(variant)

    if not is_scorable_variant_class(record.variant_class):
        return ProcessedVariantResult(
            record=record,
            score_result=None,
            status=VariantStatus.SKIPPED_UNSUPPORTED,
        )

    if (
        max_variants is not None
        and scored_count >= max_variants
    ):
        return ProcessedVariantResult(
            record=record,
            score_result=None,
            status=VariantStatus.SKIPPED_LIMIT,
        )

    try:
        score_result = compute_alphagenome_score(
            chrom=record.alphagenome_chromosome,
            pos=record.position,
            ref=record.reference,
            alt=record.alternate,
            dna_model=dna_model,
            organism=organism,
            sequence_length=sequence_length,
            selected_scorers=selected_scorers,
        )

    except Exception as exc:
        print(
            "AlphaGenome scoring failed for "
            f"{record.alphagenome_chromosome}:"
            f"{record.position} "
            f"{record.reference}>{record.alternate}: {exc}",
            file=sys.stderr,
        )
        traceback.print_exc()

        return ProcessedVariantResult(
            record=record,
            score_result=None,
            status=VariantStatus.ERROR,
        )

    return ProcessedVariantResult(
        record=record,
        score_result=score_result,
        status=VariantStatus.SCORED,
    )

def value_or_missing(value):
    """Return the configured placeholder when a value is missing."""

    return (
        value
        if value is not None
        else MISSING_VALUE
    )

def create_variant_output_values(
    score_result: AlphaGenomeScoreResult | None,
) -> VariantOutputValues:
    """Convert scoring results into VCF and TSV output values."""

    if score_result is None:
        overall = chromatin = expression = empty_score_summary()
    else:
        overall = score_result.overall_summary
        chromatin = score_result.group_summaries[CHROMATIN_GROUP.key]
        expression = score_result.group_summaries[EXPRESSION_GROUP.key]

    return VariantOutputValues(
        max_abs_score=value_or_missing(overall.max_abs_score),
        top_score=value_or_missing(overall.score),
        top_output_type=value_or_missing(overall.output_type),
        top_track=value_or_missing(overall.track),
        top_gene=value_or_missing(overall.gene),
        chromatin_score=value_or_missing(chromatin.score),
        chromatin_abs_score=value_or_missing(chromatin.max_abs_score),
        chromatin_output_type=value_or_missing(chromatin.output_type),
        chromatin_track=value_or_missing(chromatin.track),
        expression_score=value_or_missing(expression.score),
        expression_abs_score=value_or_missing(expression.max_abs_score),
        expression_gene=value_or_missing(expression.gene),
        expression_track=value_or_missing(expression.track),
    )

def annotate_variant_record(
    variant,
    *,
    output_values: VariantOutputValues,
    status: VariantStatus,
    variant_class: VariantClass,
) -> None:
    """Add AlphaGenome result fields to one VCF record."""

    values = {
        InfoField.MAX_ABS_SCORE: output_values.max_abs_score,
        InfoField.TOP_SCORE: output_values.top_score,
        InfoField.STATUS: status.value,
        InfoField.VARIANT_CLASS: variant_class.value,
        InfoField.TOP_OUTPUT_TYPE: output_values.top_output_type,
        InfoField.TOP_TRACK: output_values.top_track,
        InfoField.TOP_GENE: output_values.top_gene,
        InfoField.CHROMATIN_SCORE: output_values.chromatin_score,
        InfoField.CHROMATIN_ABS_SCORE: output_values.chromatin_abs_score,
        InfoField.CHROMATIN_OUTPUT_TYPE: output_values.chromatin_output_type,
        InfoField.CHROMATIN_TRACK: output_values.chromatin_track,
        InfoField.EXPRESSION_SCORE: output_values.expression_score,
        InfoField.EXPRESSION_ABS_SCORE: output_values.expression_abs_score,
        InfoField.EXPRESSION_GENE: output_values.expression_gene,
        InfoField.EXPRESSION_TRACK: output_values.expression_track,
    }

    for field, value in values.items():
        set_info_value(
            variant,
            field.value,
            None if value == MISSING_VALUE else value,
        ) 

def create_summary_row(
    record: VariantRecordData,
    output_values: VariantOutputValues,
    status: VariantStatus,
) -> dict[str, object]:
    """Create one compact per-variant summary row."""

    return {
        "variant_id": record.variant_id,
        "CHROM": record.original_chromosome,
        "POS": record.position,
        "REF": record.reference,
        "ALT": record.alternate,
        InfoField.MAX_ABS_SCORE.value:
            output_values.max_abs_score,
        InfoField.TOP_SCORE.value:
            output_values.top_score,
        InfoField.STATUS.value:
            status.value,
        InfoField.VARIANT_CLASS.value:
            record.variant_class.value,
        InfoField.TOP_OUTPUT_TYPE.value:
            output_values.top_output_type,
        InfoField.TOP_TRACK.value:
            output_values.top_track,
        InfoField.TOP_GENE.value:
            output_values.top_gene,
        InfoField.CHROMATIN_SCORE.value:
            output_values.chromatin_score,
        InfoField.CHROMATIN_ABS_SCORE.value:
            output_values.chromatin_abs_score,
        InfoField.CHROMATIN_OUTPUT_TYPE.value:
            output_values.chromatin_output_type,
        InfoField.CHROMATIN_TRACK.value:
            output_values.chromatin_track,
        InfoField.EXPRESSION_SCORE.value:
            output_values.expression_score,
        InfoField.EXPRESSION_ABS_SCORE.value:
            output_values.expression_abs_score,
        InfoField.EXPRESSION_GENE.value:
            output_values.expression_gene,
        InfoField.EXPRESSION_TRACK.value:
            output_values.expression_track,
    }

def run(
    vcf_path,
    output_vcf,
    output_summary_scores,
    output_tidy_scores,
    dna_model,
    organism,
    sequence_length,
    selected_scorers,
    max_variants,
):
    """Run AlphaGenome scoring and write the annotated output files."""

    summary_rows = []
    tidy_score_frames = []
    scored_count = 0

    with open_annotated_vcf(
        vcf_path,
        output_vcf,
    ) as (vcf_reader, vcf_writer):

        for variant in vcf_reader:
            processed_variant = process_variant_record(
                variant,
                dna_model=dna_model,
                organism=organism,
                sequence_length=sequence_length,
                selected_scorers=selected_scorers,
                scored_count=scored_count,
                max_variants=max_variants,
            )
            record = processed_variant.record
            score_result = processed_variant.score_result
            status = processed_variant.status

            output_values = create_variant_output_values(
                score_result
            )

            if (
                score_result is not None
                and not score_result.tidy_scores.empty
            ):
                tidy_score_frames.append(
                    score_result.tidy_scores
                )

            if status is VariantStatus.SCORED:
                scored_count += 1
                    

            annotate_variant_record(
                variant,
                output_values=output_values,
                status=status,
                variant_class=record.variant_class,
                
            )

            vcf_writer.write_record(variant)
            summary_rows.append(
                create_summary_row(record, output_values, status)
            )
    write_summary_scores( output_summary_scores, summary_rows, ) 
    write_tidy_scores( output_tidy_scores, tidy_score_frames, )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf", required=True)
    parser.add_argument("--api-key-file", required=False)
    parser.add_argument("--organism", required=True)
    parser.add_argument("--sequence-length", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary-scores", required=True, help="Path to the compact per-variant AlphaGenome summary TSV.")
    parser.add_argument("--max-variants", type=int, default=3)
    parser.add_argument("--max-scorers", type=int, default=3)
    parser.add_argument("--tidy-scores", required=True,  help="Path to the complete AlphaGenome tidy-score TSV.")
    parser.add_argument("--scorers", required=False, default=None)

    args = parser.parse_args()

    api_key = load_api_key(args.api_key_file)

    print("Initialisiere AlphaGenome API Client...", flush=True)
    model = dna_client.create(api_key)

    organism = resolve_organism(args.organism)
    sequence_length = resolve_sequence_length(args.sequence_length)
    selected_scorers = get_selected_scorers(max_scorers =args.max_scorers, scorer_names=args.scorers,)
   
    print(f"Using {len(selected_scorers)} AlphaGenome scorer(s)", flush=True)
    print(f"Requested scorer groups: {args.scorers}", flush=True)


    run(
        vcf_path=args.vcf,
        output_vcf=args.output,
        output_summary_scores=args.summary_scores,
        output_tidy_scores=args.tidy_scores,
        dna_model=model,
        organism=organism,
        sequence_length=sequence_length,
        selected_scorers=selected_scorers,
        max_variants=args.max_variants,
    )
    print("AlphaGenome scoring completed successfully.", flush=True)
    

if __name__ == "__main__":
    main()