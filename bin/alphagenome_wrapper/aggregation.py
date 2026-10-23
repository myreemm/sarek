"""
Aggregation utilities for AlphaGenome variant scores.

This module converts detailed AlphaGenome tidy-score tables into compact,
biologically grouped summary values.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd


SCORE_COLUMN_CANDIDATES = (
    "raw_score",
    "score",
    "variant_score",
)

CHROMATIN_OUTPUT_TYPES = frozenset(
    {
        "ATAC",
        "DNASE",
    }
)

EXPRESSION_OUTPUT_TYPES = frozenset(
    {
        "RNA_SEQ",
    }
)


@dataclass(frozen=True)
class TopScoreSummary:
    """Top AlphaGenome score and its associated metadata."""

    max_abs_score: Optional[float]
    output_type: Optional[str]
    track: Optional[str]
    gene: Optional[str]


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


def summarize_tidy_scores(
    score_dataframe: pd.DataFrame,
) -> TopScoreSummary:
    """
    Summarize the AlphaGenome tidy-score table for one variant.

    The row with the maximum absolute score across all selected scorers and
    tracks is used as the top result.
    """

    if score_dataframe is None or score_dataframe.empty:
        return TopScoreSummary(
            max_abs_score=None,
            output_type=None,
            track=None,
            gene=None,
        )

    numeric_scores = get_numeric_scores(score_dataframe)

    if numeric_scores is None:
        return TopScoreSummary(
            max_abs_score=None,
            output_type=None,
            track=None,
            gene=None,
        )

    top_index = numeric_scores.abs().idxmax()
    top_row = score_dataframe.loc[top_index]

    return TopScoreSummary(
        max_abs_score=float(abs(numeric_scores.loc[top_index])),
        output_type=normalize_optional_text(
            top_row.get("output_type")
        ),
        track=normalize_optional_text(
            top_row.get("track_name")
        ),
        gene=extract_top_gene(top_row),
    )


def summarize_group_score(
    score_dataframe: pd.DataFrame,
    output_types: frozenset[str],
) -> Optional[float]:
    """
    Return the maximum absolute score for selected output modalities.

    None is returned when no matching numeric score is available.
    """

    if (
        score_dataframe is None
        or score_dataframe.empty
        or "output_type" not in score_dataframe.columns
    ):
        return None

    normalized_output_types = frozenset(
        output_type.upper()
        for output_type in output_types
    )

    output_type_values = (
        score_dataframe["output_type"]
        .astype(str)
        .str.upper()
    )

    matching_rows = score_dataframe.loc[
        output_type_values.isin(normalized_output_types)
    ]

    if matching_rows.empty:
        return None

    numeric_scores = get_numeric_scores(matching_rows)

    if numeric_scores is None:
        return None

    return float(numeric_scores.abs().max())