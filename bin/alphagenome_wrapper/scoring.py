"""
AlphaGenome variant-scoring utilities.

This module performs AlphaGenome API requests and combines the detailed
results with the configured aggregation functions.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import pandas as pd

# AlphaGenome imports
from alphagenome.data import genome
from alphagenome.models import variant_scorers
from alphagenome.protos import dna_model_pb2

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



ALL_SCORERS = "all"

class ScoringMode(str, Enum):
    """Semantic scoring modes of AlphaGenome variant scorers."""

    DIFFERENTIAL = "differential"
    ACTIVE = "active"
    CONTACT_MAP = "contact_map"
    SPLICING = "splicing"
    POLYADENYLATION = "polyadenylation"


def get_scoring_mode(
    scorer,
) -> ScoringMode:
    """Determine the semantic scoring mode of one variant scorer."""

    base_scorer = scorer.base_variant_scorer

    if (
        base_scorer
        is variant_scorers.BaseVariantScorer.GENE_MASK_ACTIVE
    ):
        return ScoringMode.ACTIVE

    aggregation_type = getattr(
        scorer,
        "aggregation_type",
        None,
    )

    if aggregation_type in {
        variant_scorers.AggregationType.ACTIVE_MEAN,
        variant_scorers.AggregationType.ACTIVE_SUM,
    }:
        return ScoringMode.ACTIVE

    if base_scorer in {
        variant_scorers.BaseVariantScorer.CENTER_MASK,
        variant_scorers.BaseVariantScorer.GENE_MASK_LFC,
    }:
        return ScoringMode.DIFFERENTIAL

    if (
        base_scorer
        is variant_scorers.BaseVariantScorer.CONTACT_MAP
    ):
        return ScoringMode.CONTACT_MAP

    if base_scorer in {
        variant_scorers.BaseVariantScorer.GENE_MASK_SPLICING,
        variant_scorers.BaseVariantScorer.SPLICE_JUNCTION,
    }:
        return ScoringMode.SPLICING

    if (
        base_scorer
        is variant_scorers.BaseVariantScorer.PA_QTL
    ):
        return ScoringMode.POLYADENYLATION

    raise ValueError(
        "Unsupported AlphaGenome scorer configuration: "
        f"{scorer.name}"
    )

def to_proto_organism(
    organism,
) -> dna_model_pb2.Organism:
    """Convert the client organism enum to the AlphaGenome proto enum."""

    organism_map = {
        "HOMO_SAPIENS":
            dna_model_pb2.Organism.ORGANISM_HOMO_SAPIENS,
        "MUS_MUSCULUS":
            dna_model_pb2.Organism.ORGANISM_MUS_MUSCULUS,
    }

    try:
        return organism_map[organism.name]
    except (AttributeError, KeyError) as exc:
        raise ValueError(
            f"Unsupported AlphaGenome organism: {organism}"
        ) from exc

def get_available_scorers(
    organism,
) -> dict[str, object]:
    """Return recommended scorers supported for the selected organism."""

    proto_organism = to_proto_organism(organism)

    recommended_scorers = set(
        variant_scorers.get_recommended_scorers(
            proto_organism
        )
    )

    return {
        str(name).lower(): scorer
        for name, scorer
        in variant_scorers.RECOMMENDED_VARIANT_SCORERS.items()
        if scorer in recommended_scorers
    }



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

def resolve_requested_scorer_names(
    scorer_names: str | None,
    available_scorers: dict[str, object],
) -> list[str]:
    """Resolve the requested scorer names."""

    requested_names = parse_scorer_names(
        scorer_names
    )

    if not requested_names:
        return list(DEFAULT_SCORER_NAMES)

    if requested_names == [ALL_SCORERS]:
        return list(available_scorers)

    if ALL_SCORERS in requested_names:
        raise ValueError(
            "'all' cannot be combined with explicit scorer names."
        )

    return requested_names

def filter_signed_scores(
    score_dataframe: pd.DataFrame,
    selected_scorers,
) -> pd.DataFrame:
    """Return tidy scores produced by signed variant scorers."""

    signed_scorer_names = {
        scorer.name
        for scorer in selected_scorers
        if scorer.is_signed
    }

    if (
        score_dataframe.empty
        or "variant_scorer" not in score_dataframe.columns
    ):
        return score_dataframe.iloc[0:0].copy()

    return score_dataframe.loc[
        score_dataframe["variant_scorer"].isin(
            signed_scorer_names
        )
    ].copy()

def get_selected_scorers(
    organism,
    max_scorers: int | None = 3,
    scorer_names: str | None = None,
):
    """
    Select organism-specific recommended AlphaGenome scorers.

    If no scorer names are supplied, the validated default selection is
    used. The special value ``all`` selects every recommended scorer
    available for the organism.

    Duplicate names are removed. A positive ``max_scorers`` value acts as
    a validation limit; scorer selections are never silently truncated.
    """

    available_scorers = get_available_scorers(organism)


    requested_names = resolve_requested_scorer_names(
    scorer_names,
    available_scorers,
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
    summary_dataframe = filter_signed_scores(
    score_dataframe,
    selected_scorers,
    )

    
    overall_summary=summarize_tidy_scores(summary_dataframe)
    group_summaries = {
            group.key: summarize_group_scores(
                summary_dataframe,
                group.output_types,
            )
            for group in SCORE_GROUPS
        }

    return AlphaGenomeScoreResult(
        
        tidy_scores=score_dataframe,
        overall_summary=overall_summary,
        group_summaries=group_summaries,
    )


