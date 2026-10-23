"""
Variant parsing and classification utilities for AlphaGenome.

This module extracts variant data from cyvcf2 records, normalizes chromosome
names for AlphaGenome, and classifies variants using minimal allele
representations without modifying the original VCF records.
"""
# Imports:
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


VALID_DNA_BASES = frozenset("ACGTN")

MITOCHONDRIAL_ALIASES = frozenset(
    {
        "M",
        "MT",
        "chrM",
        "chrMT",
    }
)

# Enums:
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


SCORABLE_VARIANT_CLASSES = frozenset(
    {
        VariantClass.SNV,
        VariantClass.INSERTION,
        VariantClass.DELETION,
    }
)

# Dataclasses:
@dataclass(frozen=True)
class TrimmedAlleles:
    """Minimal allele representation used for variant classification."""

    reference: str
    alternate: str


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


def resolve_alphagenome_chromosome(
    chromosome: str,
) -> str:
    """
    Resolve a supported chromosome name for AlphaGenome.

    The chromosome name stored in the original VCF record is not modified.
    """

    chromosome_name = str(chromosome).strip()

    if not chromosome_name:
        raise ValueError(
            "The VCF record contains an empty chromosome name."
        )

    if chromosome_name in MITOCHONDRIAL_ALIASES:
        return "chrM"

    if chromosome_name.startswith("chr"):
        return chromosome_name

    return f"chr{chromosome_name}"


def is_multiallelic_variant(alternate: str) -> bool:
    """Return whether an alternate allele contains multiple alleles."""

    return "," in str(alternate)


def is_symbolic_variant(
    reference: str,
    alternate: str,
) -> bool:
    """Return whether either allele uses symbolic VCF notation."""

    reference = str(reference)
    alternate = str(alternate)

    return (
        reference.startswith("<")
        or reference.endswith(">")
        or alternate.startswith("<")
        or alternate.endswith(">")
    )


def is_valid_dna_allele(allele: str) -> bool:
    """Return whether an allele is non-empty and contains valid DNA bases."""

    normalized_allele = str(allele).upper()

    return (
        bool(normalized_allele)
        and set(normalized_allele).issubset(VALID_DNA_BASES)
    )


def trim_common_allele_context(
    reference: str,
    alternate: str,
) -> TrimmedAlleles:
    """
    Remove shared suffix and prefix bases from a biallelic variant.

    The resulting alleles are used only for classification. This function
    does not perform reference-aware left alignment and does not modify the
    original VCF record.
    """

    trimmed_reference = reference.upper()
    trimmed_alternate = alternate.upper()

    if trimmed_reference == trimmed_alternate:
        return TrimmedAlleles(
            reference="",
            alternate="",
        )

    while (
        trimmed_reference
        and trimmed_alternate
        and trimmed_reference[-1] == trimmed_alternate[-1]
    ):
        trimmed_reference = trimmed_reference[:-1]
        trimmed_alternate = trimmed_alternate[:-1]

    while (
        trimmed_reference
        and trimmed_alternate
        and trimmed_reference[0] == trimmed_alternate[0]
    ):
        trimmed_reference = trimmed_reference[1:]
        trimmed_alternate = trimmed_alternate[1:]

    return TrimmedAlleles(
        reference=trimmed_reference,
        alternate=trimmed_alternate,
    )


def classify_variant_for_alphagenome(
    reference: str,
    alternate: str,
) -> VariantClass:
    """
    Classify a biallelic variant using a minimal allele representation.

    Classification is independent of whether the current implementation
    supports scoring the resulting variant class.
    """

    normalized_reference = str(reference).strip().upper()
    normalized_alternate = str(alternate).strip().upper()

    if is_multiallelic_variant(normalized_alternate):
        return VariantClass.MULTIALLELIC

    if is_symbolic_variant(
        normalized_reference,
        normalized_alternate,
    ):
        return VariantClass.SYMBOLIC

    if (
        not is_valid_dna_allele(normalized_reference)
        or not is_valid_dna_allele(normalized_alternate)
    ):
        return VariantClass.SYMBOLIC

    if normalized_reference == normalized_alternate:
        return VariantClass.NO_CHANGE

    trimmed_alleles = trim_common_allele_context(
        normalized_reference,
        normalized_alternate,
    )

    trimmed_reference = trimmed_alleles.reference
    trimmed_alternate = trimmed_alleles.alternate

    if not trimmed_reference and not trimmed_alternate:
        return VariantClass.NO_CHANGE

    if not trimmed_reference:
        return VariantClass.INSERTION

    if not trimmed_alternate:
        return VariantClass.DELETION

    if (
        len(trimmed_reference) == 1
        and len(trimmed_alternate) == 1
    ):
        return VariantClass.SNV

    if len(trimmed_reference) == len(trimmed_alternate):
        return VariantClass.MNP

    return VariantClass.COMPLEX_DELINS


def is_scorable_variant_class(
    variant_class: VariantClass,
) -> bool:
    """Return whether the implementation currently scores a variant class."""

    return variant_class in SCORABLE_VARIANT_CLASSES


def extract_variant_record_data(
    variant,
) -> VariantRecordData:
    """Extract and classify relevant values from one cyvcf2 record."""

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