#!/usr/bin/env python3

import argparse
import csv
import os
import sys
import pandas as pd
import traceback

## AlphaGenome wrapper imports

from alphagenome_wrapper.result_io import (
    write_summary_scores,
    write_tidy_scores,
)

from alphagenome_wrapper.vcf_io import (
    INFO_CHROMATIN_SCORE,
    INFO_MAX_ABS_SCORE,
    INFO_STATUS,
    INFO_TOP_GENE,
    INFO_TOP_OUTPUT_TYPE,
    INFO_TOP_TRACK,
    INFO_VARIANT_CLASS,
    INFO_EXPRESSION_SCORE,
    open_annotated_vcf,
    set_info_value,
)
from alphagenome_wrapper.aggregation import (
    CHROMATIN_OUTPUT_TYPES,
    summarize_group_score,
    summarize_tidy_scores,
)

from alphagenome_wrapper.scoring import (
    AlphaGenomeScoreResult,
    compute_alphagenome_score,
    get_selected_scorers,
)


from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional


# AlphaGenome specific imports
from alphagenome.models import dna_client



# Constants:

STATUS_SCORED = "SCORED"
STATUS_ERROR = "ERROR"
STATUS_SKIPPED_UNSUPPORTED = "SKIPPED_UNSUPPORTED_VARIANT"
STATUS_SKIPPED_LIMIT = "SKIPPED_MAX_VARIANTS"

# Missing value used in VCF and TSV outputs
MISSING_VALUE = "."

# Enums definitions:

class VariantStatus(str, Enum):
    """Processing states assigned by the AlphaGenome wrapper."""

    SCORED = "SCORED"
    ERROR = "ERROR"
    SKIPPED_UNSUPPORTED = "SKIPPED_UNSUPPORTED_VARIANT"
    SKIPPED_LIMIT = "SKIPPED_MAX_VARIANTS"

class VariantClass(str, Enum):
    """Variant classes recognized by the AlphaGenome wrapper."""

    SNV = "SNV"
    INSERTION = "INSERTION"
    DELETION = "DELETION"
    MNP = "MNP"
    COMPLEX_DELINS = "COMPLEX_DELINS"
    MULTIALLELIC = "MULTIALLELIC"
    SYMBOLIC = "SYMBOLIC"
    NO_CHANGE = "NO_CHANGE"

# Dataclass:


@dataclass(frozen=True)
class TrimmedAlleles:
    """Minimal allele representation used for variant classification."""

    reference: str
    alternate: str
    prefix_length: int

SCORABLE_VARIANT_CLASSES = frozenset(
    {
        VariantClass.SNV,
        VariantClass.INSERTION,
        VariantClass.DELETION,
    }
)



@dataclass(frozen=True)
class VariantOutputValues:
    """Values written to the annotated VCF and summary TSV."""

    max_abs_score: float | str
    top_output_type: str
    top_track: str
    top_gene: str
    chromatin_score: float | str
    expression_score: float | str

@dataclass(frozen=True)
class VariantRecordData:
    """Structured values extracted from one VCF record."""

    original_chromosome: str
    alphagenome_chromosome: str
    position: int
    reference: str
    alternate: str
    variant_id: str
    variant_class: VariantClass

@dataclass(frozen=True)
class ProcessedVariantResult:
    """Result of processing one VCF variant record."""

    record: VariantRecordData
    score_result: AlphaGenomeScoreResult | None
    status: VariantStatus


def resolve_alphagenome_chromosome(chromosome: str) -> str:
    """
    Resolve a supported VCF chromosome name for AlphaGenome.

    The function does not modify the chromosome name stored in the output
    VCF record.
    """

    chromosome_name = str(chromosome).strip()

    if not chromosome_name:
        raise ValueError(
            "The VCF record contains an empty chromosome name."
        )

    mitochondrial_aliases = {
        "M",
        "MT",
        "chrM",
        "chrMT",
    }

    if chromosome_name in mitochondrial_aliases:
        return "chrM"

    if chromosome_name.startswith("chr"):
        return chromosome_name

    return f"chr{chromosome_name}"

def is_multiallelic_variant(alt):
    return "," in str(alt)

def is_symbolic_variant(ref, alt):
    ref = str(ref)
    alt = str(alt)
    return( 
        ref.startswith("<") 
        or ref.endswith(">") 
        or alt.startswith("<") 
        or alt.endswith(">")
    )



def is_snv(ref, alt):
    return len(ref) == 1 and len(alt) == 1


def is_simple_indel(ref, alt):
    """
    simple indel is defined as a single base insertion or deletion, where one allele is empty and the other is a single base.
    - no symbolic alleles
    - only A/C/G/T/N symbols

    """
    allowed = set("ACGTNacgtn")

    if not set(ref).issubset(allowed):
        return False

    if not set(alt).issubset(allowed):
        return False

    return len(ref) != len(alt)

def is_valid_dna_allele(allele):
    allowed = set("ACGTNacgtn")
    return bool(allele) and set(str(allele)).issubset(allowed)

def trim_common_allele_context(
    ref: str,
    alt: str,
) -> TrimmedAlleles:
    """
    Remove shared suffix and prefix bases from a biallelic variant.

    The function returns a minimal allele representation for classification
    only. It does not perform reference-aware left alignment and does not
    modify the original VCF record.

    Args:
        ref:
            Original reference allele.
        alt:
            Original alternate allele.

    Returns:
        TrimmedAlleles containing the minimal reference and alternate
        sequences and the number of removed leading bases.
    """

    trimmed_ref = ref.upper()
    trimmed_alt = alt.upper()

    if trimmed_ref == trimmed_alt:
        return TrimmedAlleles(
            reference="",
            alternate="",
            prefix_length=0,
        )

    # Remove the shared suffix first.
    while (
        trimmed_ref
        and trimmed_alt
        and trimmed_ref[-1] == trimmed_alt[-1]
    ):
        trimmed_ref = trimmed_ref[:-1]
        trimmed_alt = trimmed_alt[:-1]

    prefix_length = 0

    # Remove the shared prefix.
    while (
        trimmed_ref
        and trimmed_alt
        and trimmed_ref[0] == trimmed_alt[0]
    ):
        trimmed_ref = trimmed_ref[1:]
        trimmed_alt = trimmed_alt[1:]
        prefix_length += 1

    return TrimmedAlleles(
        reference=trimmed_ref,
        alternate=trimmed_alt,
        prefix_length=prefix_length,
    )

def classify_variant_for_alphagenome(ref: str, alt: str,) -> VariantClass:
    """
    Classify a biallelic VCF variant using a trimmed allele representation.

    Classification is independent of whether the current prototype scores
    the resulting variant class.
    """

    normalized_ref = str(ref).strip().upper()
    normalized_alt = str(alt).strip().upper()

    if is_multiallelic_variant(normalized_alt):
        return VariantClass.MULTIALLELIC

    if is_symbolic_variant(normalized_ref, normalized_alt):
        return VariantClass.SYMBOLIC

    if (
        not is_valid_dna_allele(normalized_ref)
        or not is_valid_dna_allele(normalized_alt)
    ):
        return VariantClass.SYMBOLIC

    if normalized_ref == normalized_alt:
        return VariantClass.NO_CHANGE

    trimmed = trim_common_allele_context(
        normalized_ref,
        normalized_alt,
    )

    trimmed_ref = trimmed.reference
    trimmed_alt = trimmed.alternate

    if not trimmed_ref and not trimmed_alt:
        return VariantClass.NO_CHANGE

    if not trimmed_ref:
        return VariantClass.INSERTION

    if not trimmed_alt:
        return VariantClass.DELETION

    if len(trimmed_ref) == 1 and len(trimmed_alt) == 1:
        return VariantClass.SNV

    if len(trimmed_ref) == len(trimmed_alt):
        return VariantClass.MNP

    return VariantClass.COMPLEX_DELINS
def is_scorable_variant_class(variant_class: VariantClass,) -> bool:
    """Return whether a variant class is currently scored by the prototype."""

    return variant_class in SCORABLE_VARIANT_CLASSES

 
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


def sanitize_info_value(value):
    """
    Make a value safe for VCF-INFO fields.
    Removes characters that could make INFO fields illegible or malformed.

    """

    if value is None:
        return "."

    value = str(value)

    if value == "" or value.lower() == "nan":
        return "."

    return (
        value
        .replace(";", "_")
        .replace(",", "_")
        .replace("\t", "_")
        .replace(" ", "_")
    )

def extract_variant_record_data(variant) -> VariantRecordData:
    """
    Extract and classify the relevant values from one cyvcf2 record.

    The original chromosome name is preserved for the output VCF. A separate
    chromosome name is resolved for the AlphaGenome API request.
    """

    original_chromosome = str(variant.CHROM).strip()
    alphagenome_chromosome = resolve_alphagenome_chromosome(
        original_chromosome
    )

    position = int(variant.POS)
    reference = str(variant.REF).upper()

    alternate_alleles = [
        str(allele).upper()
        for allele in variant.ALT
        if allele is not None
    ]

    if not alternate_alleles:
        alternate = ""
        variant_class = VariantClass.NO_CHANGE

    elif len(alternate_alleles) > 1:
        alternate = ",".join(alternate_alleles)
        variant_class = VariantClass.MULTIALLELIC

    else:
        alternate = alternate_alleles[0]
        variant_class = classify_variant_for_alphagenome(
            reference,
            alternate,
        )

    variant_id = (
        f"{original_chromosome}_{position}_{reference}_{alternate}"
    )

    return VariantRecordData(
        original_chromosome=original_chromosome,
        alphagenome_chromosome=alphagenome_chromosome,
        position=position,
        reference=reference,
        alternate=alternate,
        variant_id=variant_id,
        variant_class=variant_class,
    )

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

def create_variant_output_values(
    score_result: AlphaGenomeScoreResult | None,
) -> VariantOutputValues:
    """
    Convert an optional score result into values for VCF and TSV output.
    """

    if score_result is None:
        return VariantOutputValues(
            max_abs_score=MISSING_VALUE,
            top_output_type=MISSING_VALUE,
            top_track=MISSING_VALUE,
            top_gene=MISSING_VALUE,
            chromatin_score=MISSING_VALUE,
            expression_score=MISSING_VALUE,
        )

    return VariantOutputValues(
        max_abs_score=(
            score_result.max_abs_score
            if score_result.max_abs_score is not None
            else MISSING_VALUE
        ),
        top_output_type=(
            score_result.top_output_type or MISSING_VALUE
        ),
        top_track=(
            score_result.top_track or MISSING_VALUE
        ),
        top_gene=(
            score_result.top_gene or MISSING_VALUE
        ),
        chromatin_score=(
            score_result.chromatin_score
            if score_result.chromatin_score is not None
            else MISSING_VALUE
        ),
        expression_score=(
            score_result.expression_score
            if score_result.expression_score is not None
            else MISSING_VALUE
        ),
    )

def annotate_variant_record(
    variant,
    *,
    max_abs_score,
    status: VariantStatus,
    variant_class: VariantClass,
    top_output_type,
    top_track,
    top_gene,
    chromatin_score,
    expression_score,
) -> None:
    """
    Add AlphaGenome result fields to a cyvcf2 variant record.

    Missing optional values are omitted from the INFO column. Processing
    status and variant class are always written.
    """

    set_info_value(
        variant,
        INFO_MAX_ABS_SCORE,
        None
        if max_abs_score == MISSING_VALUE
        else max_abs_score,
    )
    set_info_value(
        variant,
        INFO_STATUS,
        status.value,
    )
    set_info_value(
        variant,
        INFO_VARIANT_CLASS,
        variant_class.value,
    )
    set_info_value(
        variant,
        INFO_TOP_OUTPUT_TYPE,
        None
        if top_output_type == MISSING_VALUE
        else str(top_output_type),
    )
    set_info_value(
        variant,
        INFO_TOP_TRACK,
        None
        if top_track == MISSING_VALUE
        else str(top_track),
    )
    set_info_value(
        variant,
        INFO_TOP_GENE,
        None
        if top_gene == MISSING_VALUE
        else str(top_gene),
    )
    set_info_value(
        variant,
        INFO_CHROMATIN_SCORE,
        None
        if chromatin_score == MISSING_VALUE
        else chromatin_score,
    )
    set_info_value(
        variant,
        INFO_EXPRESSION_SCORE,
        None
        if expression_score == MISSING_VALUE
        else expression_score,
        )



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
                max_abs_score=output_values.max_abs_score,
                status=status,
                variant_class=record.variant_class,
                top_output_type=output_values.top_output_type,
                top_track=output_values.top_track,
                top_gene=output_values.top_gene,
                chromatin_score=output_values.chromatin_score,
                expression_score=output_values.expression_score,
            )

            vcf_writer.write_record(variant)
            summary_rows.append({"variant_id": record.variant_id,"CHROM": record.original_chromosome,"POS": record.position,"REF": record.reference,"ALT": record.alternate,
                INFO_MAX_ABS_SCORE:( output_values.max_abs_score),
                INFO_STATUS: status.value,
                INFO_VARIANT_CLASS: (
                        record.variant_class.value),
                INFO_TOP_OUTPUT_TYPE: (
                        output_values.top_output_type),
                INFO_TOP_TRACK: (
                        output_values.top_track),
                INFO_TOP_GENE: (
                        output_values.top_gene),
                INFO_CHROMATIN_SCORE: (
                        output_values.chromatin_score),
                INFO_EXPRESSION_SCORE: (
                        output_values.expression_score),
            }
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
    parser.add_argument("--max-scorers", type=int, default=1)
    parser.add_argument("--tidy-scores", required=True,  help="Path to the complete AlphaGenome tidy-score TSV.")
    parser.add_argument("--scorers", required=False, default=None)

    args = parser.parse_args()

    api_key = load_api_key(args.api_key_file)

    print("Initialisiere AlphaGenome API Client...", flush=True)
    model = dna_client.create(api_key)

    organism = resolve_organism(args.organism)
    sequence_length = resolve_sequence_length(args.sequence_length)
    selected_scorers = get_selected_scorers(max_scorers =args.max_scorers, scorer_names=args.scorers,)
    #debug print to check the number of selected scorers
    print(
        f"Requested scorer keys: {args.scorers}",
        file=sys.stderr,
    )

    print(
        f"Selected scorer count: {len(selected_scorers)}",
        file=sys.stderr,
    )

    for scorer in selected_scorers:
        print(
            f"Selected scorer: {scorer.name}",
            file=sys.stderr,
        )


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
    print("Initializing AlphaGenome API Client...")
    print(f"Using {len(selected_scorers)} AlphaGenome scorer(s)", flush=True) # temporary debug print
    

if __name__ == "__main__":
    main()