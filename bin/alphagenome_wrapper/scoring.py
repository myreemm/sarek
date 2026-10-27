"""
AlphaGenome variant-scoring utilities.

This module performs AlphaGenome API requests and combines the detailed
results with the configured aggregation functions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

# AlphaGenome imports
from alphagenome.data import genome
from alphagenome.models import variant_scorers

# AlphaGenome wrapper imports

from alphagenome_wrapper.aggregation import (
    SCORE_GROUPS,
    TopScoreSummary,
    summarize_group_scores,
    summarize_tidy_scores,
)
# Constants:

DEFAULT_SCORER_NAMES = (
    "atac",
    "dnase",
    "rna_seq",
)

ACTIVE_SCORER_SUFFIX = "_active"



@dataclass(frozen=True)
class AlphaGenomeScoreResult:
    """
    Aggregated result of an AlphaGenome variant-scoring request.

    Attributes:
        tidy_scores:
            Complete tidy AlphaGenome results for the scored variant.
        overall_summary:
            Signed top score and associated context across all selected
            scorer results.
        chromatin_summary:
            Signed top score and associated context across ATAC and
            DNASE results.
        expression_summary:
            Signed top score and associated context across RNA_SEQ
            results.
    """

    tidy_scores: pd.DataFrame
    overall_summary: TopScoreSummary
    group_summaries: dict[str, TopScoreSummary]
    



def parse_scorer_names(
    scorer_names: str | None,
) -> list[str]:
    """
    Parse comma-separated AlphaGenome scorer names.

    Duplicate names are removed while preserving their requested order.
    """

    normalized_text = str(
        scorer_names or ""
    ).strip().lower()

    if normalized_text in {
        "",
        "none",
        "null",
    }:
        return []

    parsed_names = [
        name.strip()
        for name in normalized_text.split(",")
        if name.strip()
    ]

    return list(
        dict.fromkeys(parsed_names)
    )


def get_selected_scorers(
    max_scorers: int | None = 3,
    scorer_names: str | None = None,
):
    """
    Select validated differential AlphaGenome scorer configurations.

    If no explicit scorer names are supplied, the wrapper uses ATAC,
    DNASE, and RNA_SEQ as its documented default selection.

    Duplicate scorer names are removed. A positive max_scorers value acts
    only as a validation limit. The selection is never silently truncated.

    Active-allele scorers are currently rejected because their absolute
    activity values are not compatible with the prototype's differential
    compact summaries.
    """

    available_scorers = {
        str(key).lower(): scorer
        for key, scorer
        in variant_scorers.RECOMMENDED_VARIANT_SCORERS.items()
    }

    requested_names = (
        parse_scorer_names(scorer_names)
        or list(DEFAULT_SCORER_NAMES)
    )

    unknown_names = [
        name
        for name in requested_names
        if name not in available_scorers
    ]

    if unknown_names:
        raise ValueError(
            "Unknown AlphaGenome scorer(s): "
            + ", ".join(unknown_names)
            + ". Available  scorer keys are: "
            + ", ".join(sorted(available_scorers))
        )

    active_names = [
        name
        for name in requested_names
        if name.endswith(ACTIVE_SCORER_SUFFIX)
    ]

    if active_names:
        raise ValueError(
            "Active-allele scorers are not yet supported by "
            "the compact differential summary: "
            + ", ".join(active_names)
        )

    if (
        max_scorers is not None
        and max_scorers > 0
        and len(requested_names) > max_scorers
    ):
        raise ValueError(
            "The number of selected AlphaGenome scorers "
            f"({len(requested_names)}) exceeds the configured "
            f"maximum ({max_scorers}). Selected scorers: "
            + ", ".join(requested_names)
        )

    return [
        available_scorers[name]
        for name in requested_names
    ]


def compute_alphagenome_score(chrom, pos, ref, alt, dna_model, organism, sequence_length, selected_scorers,) -> AlphaGenomeScoreResult:
    """
    Score one variant and summarize its AlphaGenome results.

    Top rows are selected by maximum absolute raw score. The original
    signed score and the metadata from the selected row are preserved.
    """
    variant_id = f"{chrom}_{pos}_{ref}_{alt}"

    variant = genome.Variant(
        chromosome=str(chrom),
        position=int(pos),
        reference_bases=str(ref),
        alternate_bases=str(alt),
        name=variant_id,
    )

    interval = variant.reference_interval.resize(sequence_length)

    variant_scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=selected_scorers,
        organism=organism,
    )

    score_dataframe = variant_scorers.tidy_scores([variant_scores])
    
    overall_summary=summarize_tidy_scores(score_dataframe)
    group_summaries = {
            group.key: summarize_group_scores(
                score_dataframe,
                group.output_types,
            )
            for group in SCORE_GROUPS
        }

    return AlphaGenomeScoreResult(
        
        tidy_scores=score_dataframe,
        overall_summary=overall_summary,
        group_summaries=group_summaries,
    )


