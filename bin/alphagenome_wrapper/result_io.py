"""
Output utilities for AlphaGenome result tables.

This module centralizes the column definitions and writing logic for the
compact summary table and the complete tidy-score table.
"""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from alphagenome_wrapper.vcf_io import InfoField




SUMMARY_COLUMNS = (
    "variant_id",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    InfoField.TOP_SCORE.value,
    InfoField.MAX_ABS_SCORE.value,
    InfoField.STATUS.value,
    InfoField.VARIANT_CLASS.value,
    InfoField.TOP_OUTPUT_TYPE.value,
    InfoField.TOP_TRACK.value,
    InfoField.TOP_GENE.value,
    InfoField.CHROMATIN_SCORE.value,
    InfoField.CHROMATIN_ABS_SCORE.value,
    InfoField.CHROMATIN_OUTPUT_TYPE.value,
    InfoField.CHROMATIN_TRACK.value,
    InfoField.EXPRESSION_SCORE.value,
    InfoField.EXPRESSION_ABS_SCORE.value,
    InfoField.EXPRESSION_GENE.value,
    InfoField.EXPRESSION_TRACK.value,
)


EMPTY_TIDY_COLUMNS = (
    "variant_id",
    "output_type",
    "variant_scorer",
    "track_name",
    "raw_score",
)


def write_summary_scores(
    output_path: str | Path,
    summary_rows: Sequence[Mapping[str, Any]],
) -> None:
    """
    Write the compact per-variant AlphaGenome summary table.

    Args:
        output_path:
            Destination path for the tab-separated summary file.
        summary_rows:
            Per-variant mappings whose keys correspond to SUMMARY_COLUMNS.
    """

    output_path = Path(output_path)

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as output_handle:
        writer = csv.DictWriter(
            output_handle,
            fieldnames=SUMMARY_COLUMNS,
            delimiter="\t",
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)


def write_tidy_scores(
    output_path: str | Path,
    tidy_score_frames: Sequence[pd.DataFrame],
) -> None:
    """
    Write the complete AlphaGenome tidy-score table.

    An empty table with a stable minimum schema is written when no variant
    produced tidy-score results.

    Args:
        output_path:
            Destination path for the tab-separated tidy-score file.
        tidy_score_frames:
            AlphaGenome tidy-score tables collected during the run.
    """

    output_path = Path(output_path)

    non_empty_frames = [
        dataframe
        for dataframe in tidy_score_frames
        if dataframe is not None and not dataframe.empty
    ]

    if non_empty_frames:
        tidy_dataframe = pd.concat(
            non_empty_frames,
            ignore_index=True,
        )
    else:
        tidy_dataframe = pd.DataFrame(
            columns=EMPTY_TIDY_COLUMNS
        )

    tidy_dataframe.to_csv(
        output_path,
        sep="\t",
        index=False,
    )
