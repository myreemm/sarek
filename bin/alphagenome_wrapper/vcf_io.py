"""
VCF input and output utilities for the AlphaGenome wrapper.

This module centralizes VCF header definitions and provides helpers for
opening, annotating, and writing VCF records with cyvcf2.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import closing, contextmanager, suppress
from pathlib import Path
from typing import Any

from cyvcf2 import VCF, Writer
from enum import StrEnum


# Grouped AlphaGenome INFO fields

class InfoField(StrEnum):
    """VCF INFO field identifiers emitted by the AlphaGenome wrapper."""

    MAX_ABS_SCORE = "AG_MAX_ABS_SCORE"
    TOP_SCORE = "AG_TOP_SCORE"

    STATUS = "AG_STATUS"
    VARIANT_CLASS = "AG_VARIANT_CLASS"

    TOP_OUTPUT_TYPE = "AG_TOP_OUTPUT_TYPE"
    TOP_TRACK = "AG_TOP_TRACK"
    TOP_GENE = "AG_TOP_GENE"

    CHROMATIN_SCORE = "AG_CHROMATIN_SCORE"
    CHROMATIN_ABS_SCORE = "AG_CHROMATIN_ABS_SCORE"
    CHROMATIN_OUTPUT_TYPE = "AG_CHROMATIN_OUTPUT_TYPE"
    CHROMATIN_TRACK = "AG_CHROMATIN_TRACK"

    EXPRESSION_SCORE = "AG_EXPRESSION_SCORE"
    EXPRESSION_ABS_SCORE = "AG_EXPRESSION_ABS_SCORE"
    EXPRESSION_GENE = "AG_EXPRESSION_GENE"
    EXPRESSION_TRACK = "AG_EXPRESSION_TRACK"

# One row per INFO field: (ID, VCF type, description).
# All fields are single-valued, so Number is "1" throughout.
_INFO_FIELD_TABLE: tuple[tuple[InfoField, str, str], ...] = (
    (
        InfoField.MAX_ABS_SCORE,
        "Float",
        "Absolute magnitude of the global top differential AlphaGenome score",
    ),
    (
        InfoField.TOP_SCORE,
        "Float",
        "Signed global top differential AlphaGenome score",
    ),
    (InfoField.STATUS, "String", "AlphaGenome processing status"),
    (
        InfoField.VARIANT_CLASS,
        "String",
        "Variant class assigned by the AlphaGenome wrapper",
    ),
    (
        InfoField.TOP_OUTPUT_TYPE,
        "String",
        "Output type associated with AG_TOP_SCORE",
    ),
    (InfoField.TOP_TRACK, "String", "Track associated with AG_TOP_SCORE"),
    (
        InfoField.TOP_GENE,
        "String",
        "Gene associated with AG_TOP_SCORE, if available",
    ),
    (
        InfoField.CHROMATIN_SCORE,
        "Float",
        "Signed top chromatin-accessibility score selected "
        "from ATAC and DNASE results",
    ),
    (
        InfoField.CHROMATIN_ABS_SCORE,
        "Float",
        "Absolute magnitude of AG_CHROMATIN_SCORE",
    ),
    (
        InfoField.CHROMATIN_OUTPUT_TYPE,
        "String",
        "ATAC or DNASE output type associated with AG_CHROMATIN_SCORE",
    ),
    (
        InfoField.CHROMATIN_TRACK,
        "String",
        "Track associated with AG_CHROMATIN_SCORE",
    ),
    (
        InfoField.EXPRESSION_SCORE,
        "Float",
        "Signed top RNA_SEQ expression score",
    ),
    (
        InfoField.EXPRESSION_ABS_SCORE,
        "Float",
        "Absolute magnitude of AG_EXPRESSION_SCORE",
    ),
    (
        InfoField.EXPRESSION_GENE,
        "String",
        "Gene associated with AG_EXPRESSION_SCORE",
    ),
    (
        InfoField.EXPRESSION_TRACK,
        "String",
        "Track associated with AG_EXPRESSION_SCORE",
    ),
)

# Header definitions in the dict format expected by cyvcf2.
ALPHAGENOME_INFO_FIELDS: tuple[dict[str, str], ...] = tuple(
    {
        "ID": field_id.value,
        "Description": description,
        "Type": field_type,
        "Number": "1",
    }
    for field_id, field_type, description in _INFO_FIELD_TABLE
)


def add_alphagenome_info_headers(vcf_reader: VCF) -> None:
    """
    Register all AlphaGenome INFO fields in a VCF header.

    The header must be modified before the cyvcf2 Writer is created.
    Fields that already exist in the header are left untouched.
    """

    for field_definition in ALPHAGENOME_INFO_FIELDS:
        if not vcf_reader.contains(field_definition["ID"]):
            vcf_reader.add_info_to_header(field_definition)


@contextmanager
def open_annotated_vcf(
    input_path: str | Path,
    output_path: str | Path,
) -> Iterator[tuple[VCF, Writer]]:
    """
    Open an input VCF and an output writer with AlphaGenome INFO headers.

    Reader and writer are closed automatically, also when an exception
    occurs.

    Args:
        input_path: Path to an input VCF, compressed VCF, or BCF.
        output_path: Path to the annotated output VCF or BCF.

    Yields:
        Tuple containing the configured VCF reader and writer.
    """

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(f"Input VCF does not exist: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with closing(VCF(str(input_path))) as vcf_reader:
        add_alphagenome_info_headers(vcf_reader)

        with closing(Writer(str(output_path), vcf_reader)) as vcf_writer:
            yield vcf_reader, vcf_writer


def set_info_value(variant: Any, field_name: str, value: Any) -> None:
    """
    Set or clear one INFO value on a cyvcf2 variant record.

    Args:
        variant: cyvcf2 Variant record.
        field_name: Registered INFO field identifier.
        value: Value to write. None removes an existing field value.
    """

    if value is None:
        with suppress(KeyError):
            del variant.INFO[field_name]
    else:
        variant.INFO[field_name] = value