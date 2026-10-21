"""
AlphaGenome variant-scoring utilities.

This module performs AlphaGenome API requests and combines the detailed
results with the configured aggregation functions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd
import sys

# AlpgaGenome imports
from alphagenome.data import genome
from alphagenome.models import variant_scorers

from alphagenome_wrapper.aggregation import (
    CHROMATIN_OUTPUT_TYPES,
    EXPRESSION_OUTPUT_TYPES,
    summarize_group_score,
    summarize_tidy_scores,
)


@dataclass(frozen=True)
class AlphaGenomeScoreResult:

    """
    Aggregated result of an AlphaGenome variant-scoring request.

    Attributes:
        max_abs_score:
            Maximum absolute raw score across the selected scorers and tracks.
        tidy_scores:
            Complete tidy AlphaGenome results for the scored variant.
        top_output_type:
            Output modality associated with the maximum absolute score.
        top_track:
            Track associated with the maximum absolute score.
        top_gene:
            Gene associated with the maximum absolute score, if available.
        chromatin_score:
            Maximum absolute score across chromatin-related output types.
    """

    max_abs_score: Optional[float]
    tidy_scores: pd.DataFrame
    top_output_type: Optional[str]
    top_track: Optional[str]
    top_gene: Optional[str]
    chromatin_score: Optional[float]
    expression_score: Optional[float]
    


# Helper function to select the requested scorers from the available scorers
def get_selected_scorers(
    max_scorers: int = 1,
    scorer_names: str | None = None,
):
    """
    Select AlphaGenome variant scorers.

    Explicit scorer names take precedence over max_scorers. The limit is
    applied only when no explicit scorer names are provided.
    """

    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS

    normalized_scorer_names = (
        str(scorer_names).strip().lower()
        if scorer_names is not None
        else ""
    )

    if normalized_scorer_names not in {
        "",
        "none",
        "null",
    }:
        requested_names = [
            name.strip().lower()
            for name in normalized_scorer_names.split(",")
            if name.strip()
        ]

        available_scorers = {
            str(key).lower(): scorer
            for key, scorer in all_scorers.items()
        }

        unknown_names = [
            name
            for name in requested_names
            if name not in available_scorers
        ]

        if unknown_names:
            raise ValueError(
                "Unknown AlphaGenome scorer(s): "
                + ", ".join(unknown_names)
                + ". Available scorer keys are: "
                + ", ".join(sorted(available_scorers))
            )

        return [
            available_scorers[name]
            for name in requested_names
        ]

    scorers = list(all_scorers.values())

    if max_scorers is None or max_scorers <= 0:
        return scorers

    return scorers[:max_scorers]


def compute_alphagenome_score(chrom, pos, ref, alt, dna_model, organism, sequence_length, selected_scorers,) -> AlphaGenomeScoreResult:
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
    for scorer, scorer_result in zip(
        selected_scorers,
        variant_scores,
        strict=True,
    ):
        print(
            "Scorer result:",
            scorer.name,
            "shape=",
            scorer_result.shape,
        )

    score_dataframe = variant_scorers.tidy_scores([variant_scores])

    top_summary = summarize_tidy_scores(score_dataframe)

    chromatin_score = summarize_group_score(score_dataframe, CHROMATIN_OUTPUT_TYPES,)

    expression_score = summarize_group_score(score_dataframe, EXPRESSION_OUTPUT_TYPES,)

    return AlphaGenomeScoreResult(
        max_abs_score=top_summary.max_abs_score,
        tidy_scores=score_dataframe,
        top_output_type=top_summary.output_type,
        top_track=top_summary.track,
        top_gene=top_summary.gene,
        chromatin_score=chromatin_score,
        expression_score=expression_score,
    )

