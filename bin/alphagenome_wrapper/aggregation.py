"""
Aggregation utilities for AlphaGenome variant scores.

This module converts detailed AlphaGenome tidy-score tables into compact,
biologically grouped summary values.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Optional
import pandas as pd


SCORE_COLUMN_CANDIDATES = (
    "raw_score",
    "score",
    "variant_score",
)


# Dataclasses:



@dataclass(frozen=True)
class TopScoreSummary:
    """Top AlphaGenome score and its associated metadata."""

    score: Optional[float]
    max_abs_score: Optional[float]
    output_type: Optional[str]
    track: Optional[str]
    gene: Optional[str]

class ScoreFieldSuffix(StrEnum):
    """Suffixes used for compact grouped score fields."""

    SCORE = "SCORE"
    ABS_SCORE = "ABS_SCORE"
    OUTPUT_TYPE = "OUTPUT_TYPE"
    TRACK = "TRACK"
    GENE = "GENE"


@dataclass(frozen=True)
class ScoreGroup:
    """Configuration for one biological score group."""

    key: str
    output_types: frozenset[str]
    context_fields: tuple[ScoreFieldSuffix, ...]


CHROMATIN_GROUP = ScoreGroup(
    key="CHROMATIN",
    output_types=frozenset(
        {
            "ATAC",
            "DNASE",
        }
    ),
    context_fields=(
        ScoreFieldSuffix.OUTPUT_TYPE,
        ScoreFieldSuffix.TRACK,
    ),
)


EXPRESSION_GROUP = ScoreGroup(
    key="EXPRESSION",
    output_types=frozenset(
        {
            "RNA_SEQ",
        }
    ),
    context_fields=(
        ScoreFieldSuffix.TRACK,
        ScoreFieldSuffix.GENE,
    ),
)


SCORE_GROUPS = (
    CHROMATIN_GROUP,
    EXPRESSION_GROUP,
)

# Helper functions:

def empty_score_summary() -> TopScoreSummary:
    """Return an empty score summary."""

    return TopScoreSummary(
        score=None,
        max_abs_score=None,
        output_type=None,
        track=None,
        gene=None,
    )

def find_score_column(
    score_dataframe: pd.DataFrame,
) -> Optional[str]:
    """Return the first supported score column in a tidy-score table."""

    return next(
        (
            column
            for column in SCORE_COLUMN_CANDIDATES
            if column in score_dataframe.columns
        ),
        None,
    )


def get_numeric_scores(
    score_dataframe: pd.DataFrame,
) -> Optional[pd.Series]:
    """
    Return the valid numeric score values from a tidy-score table.

    None is returned when no supported score column or no numeric values
    are available.
    """

    score_column = find_score_column(score_dataframe)

    if score_column is None:
        return None

    numeric_scores = pd.to_numeric(
        score_dataframe[score_column],
        errors="coerce",
    ).dropna()

    if numeric_scores.empty:
        return None

    return numeric_scores


def normalize_optional_text(value) -> Optional[str]:
    """Return a stripped string or None for a missing or empty value."""

    if value is None or pd.isna(value):
        return None

    normalized_value = str(value).strip()

    return normalized_value or None


def extract_top_gene(
    top_row: pd.Series,
) -> Optional[str]:
    """
    Extract the most informative gene identifier from a tidy-score row.

    The gene name is preferred. The gene identifier is used as fallback.
    """

    return (
        normalize_optional_text(top_row.get("gene_name"))
        or normalize_optional_text(top_row.get("gene_id"))
    )


def summarize_scores(
    score_dataframe: pd.DataFrame,
    output_types: frozenset[str] | None = None,
) -> TopScoreSummary:
    """
    Select the row with the largest absolute score.

    The original signed score and all returned metadata are taken from
    the same tidy-score row. When output_types is provided, only matching
    output modalities are considered.
    """

    if score_dataframe is None or score_dataframe.empty:
        return empty_score_summary()

    selected_rows = score_dataframe

    if output_types is not None:
        if "output_type" not in selected_rows.columns:
            return empty_score_summary()

        normalized_output_types = frozenset(
            output_type.upper()
            for output_type in output_types
        )

        output_type_values = (
            selected_rows["output_type"]
            .astype(str)
            .str.upper()
        )

        selected_rows = selected_rows.loc[
            output_type_values.isin(
                normalized_output_types
            )
        ]

    if selected_rows.empty:
        return empty_score_summary()

    numeric_scores = get_numeric_scores(
        selected_rows
    )

    if numeric_scores is None:
        return empty_score_summary()

    top_index = numeric_scores.abs().idxmax()
    top_row = selected_rows.loc[top_index]
    signed_score = float(
        numeric_scores.loc[top_index]
    )

    return TopScoreSummary(
        score=signed_score,
        max_abs_score=abs(signed_score),
        output_type=normalize_optional_text(
            top_row.get("output_type")
        ),
        track=normalize_optional_text(
            top_row.get("track_name")
        ),
        gene=extract_top_gene(top_row),
    )

def summarize_tidy_scores(
    score_dataframe: pd.DataFrame,
) -> TopScoreSummary:
    """
    Summarize all AlphaGenome score results.

    The top row is selected by maximum absolute score while preserving
    the original signed score and the metadata from the same row.
    """

    return summarize_scores(
        score_dataframe
    )


def summarize_group_scores(
    score_dataframe: pd.DataFrame,
    output_types: frozenset[str],
) -> TopScoreSummary:
    """
    Summarize AlphaGenome scores for selected output modalities.

    The top row is selected by maximum absolute score while preserving
    the original signed score and the metadata from the same row.
    """

    return summarize_scores(
        score_dataframe,
        output_types,
    )