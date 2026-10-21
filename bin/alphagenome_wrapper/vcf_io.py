"""
VCF input and output utilities for the AlphaGenome wrapper.
This module centralizes VCF header definitions and provides helpers for
opening, annotating, and writing VCF records with cyvcf2.
"""
 
from __future__ import annotations
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from cyvcf2 import VCF, Writer



INFO_MAX_ABS_SCORE = "AG_MAX_ABS_SCORE"
INFO_STATUS = "AG_STATUS"
INFO_VARIANT_CLASS = "AG_VARIANT_CLASS"
INFO_TOP_OUTPUT_TYPE = "AG_TOP_OUTPUT_TYPE"
INFO_TOP_TRACK = "AG_TOP_TRACK"
INFO_TOP_GENE = "AG_TOP_GENE"
INFO_CHROMATIN_SCORE = "AG_CHROMATIN_SCORE"
INFO_EXPRESSION_SCORE = "AG_EXPRESSION_SCORE"


ALPHAGENOME_INFO_FIELDS: tuple[dict[str, str], ...] = (
    {
        "ID": INFO_MAX_ABS_SCORE,
        "Description": (
            "Maximum absolute AlphaGenome raw score across the selected "
            "scorers and tracks"
        ),
        "Type": "Float",
        "Number": "1",
    },
    {
        "ID": INFO_STATUS,
        "Description": "AlphaGenome processing status",
        "Type": "String",
        "Number": "1",
    },
    {
        "ID": INFO_VARIANT_CLASS,
        "Description": (
            "Variant class assigned by the AlphaGenome wrapper"
        ),
        "Type": "String",
        "Number": "1",
    },
    {
        "ID": INFO_TOP_OUTPUT_TYPE,
        "Description": (
            "Output type associated with the maximum absolute score"
        ),
        "Type": "String",
        "Number": "1",
    },
    {
        "ID": INFO_TOP_TRACK,
        "Description": (
            "Track associated with the maximum absolute score"
        ),
        "Type": "String",
        "Number": "1",
    },
    {
        "ID": INFO_TOP_GENE,
        "Description": (
            "Gene associated with the maximum absolute score, if available"
        ),
        "Type": "String",
        "Number": "1",
    },
    {
        "ID": INFO_CHROMATIN_SCORE,
        "Description": (
            "Maximum absolute chromatin-related AlphaGenome score"
        ),
        "Type": "Float",
        "Number": "1",
    },
    {
        "ID": INFO_EXPRESSION_SCORE,
        "Description": (
            "Maximum absolute expression-related AlphaGenome score"
        ),
        "Type": "Float",
        "Number": "1",
    },
)


def add_alphagenome_info_headers(vcf_reader: VCF) -> None:
    """
    Register all AlphaGenome INFO fields in a VCF header.

    The header must be modified before the cyvcf2 Writer is created.

    Args:
        vcf_reader:
            Open cyvcf2 VCF reader whose header will be extended.
    """

    for field_definition in ALPHAGENOME_INFO_FIELDS:
        field_id = field_definition["ID"]

        if not vcf_reader.contains(field_id):
            vcf_reader.add_info_to_header(field_definition)


@contextmanager
def open_annotated_vcf(
    input_path: str | Path,
    output_path: str | Path,
) -> Iterator[tuple[VCF, Writer]]:
    """
    Open an input VCF and an output writer with AlphaGenome INFO headers.

    The reader and writer are closed automatically, including when an
    exception occurs.

    Args:
        input_path:
            Path to an input VCF, compressed VCF, or BCF.
        output_path:
            Path to the annotated output VCF or BCF.

    Yields:
        Tuple containing the configured VCF reader and writer.
    """

    input_path = Path(input_path)
    output_path = Path(output_path)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Input VCF does not exist: {input_path}"
        )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    vcf_reader = VCF(str(input_path))
    vcf_writer: Writer | None = None

    try:
        add_alphagenome_info_headers(vcf_reader)
        vcf_writer = Writer(str(output_path), vcf_reader)

        yield vcf_reader, vcf_writer

    finally:
        if vcf_writer is not None:
            vcf_writer.close()

        vcf_reader.close()


def set_info_value(
    variant: Any,
    field_name: str,
    value: Any,
) -> None:
    """
    Set or clear one INFO value on a cyvcf2 variant record.

    Args:
        variant:
            cyvcf2 Variant record.
        field_name:
            Registered INFO field identifier.
        value:
            Value to write. None removes an existing field value.
    """

    if value is None:
        try:
            del variant.INFO[field_name]
        except KeyError:
            pass

        return

    variant.INFO[field_name] = value