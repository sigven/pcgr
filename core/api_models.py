from typing import Union, List, Dict, Any
from pydantic import BaseModel, Field


class SNVInput(BaseModel):
    """Model for simple nucleotide variants (SNVs)."""
    chromosome: str = Field(description="Chromosome (e.g., 'chr1', '1')")
    position: int = Field(description="Genomic position")
    reference_allele: str = Field(description="Reference allele")
    alternate_allele: str = Field(description="Alternate allele")


class CNVInput(BaseModel):
    """Model for Copy Number Variations (CNVs)."""
    gene: str = Field(description="Gene symbol")
    copy_number_change: str = Field(description="Type of copy number change (e.g., 'amplification', 'loss', 'high-level amplification')")


class SVInput(BaseModel):
    """Model for Structural Variants (SVs)."""
    sv_type: str = Field(description="Structural variant type (e.g., 'translocation', 'fusion')")
    genes: List[str] = Field(description="List of genes involved (e.g., ['BCR', 'ABL1'])")


Variant = Union[SNVInput, CNVInput, SVInput]


class ClassificationRequest(BaseModel):
    """Main model for an API classification request."""
    variant: Variant = Field(description="Variant to classify (SNV, CNV, or SV)")
    oncotree_code: str = Field(description="OncoTree disease code")


class ClassificationReport(BaseModel):
    """Final report object returned by the API."""
    request: ClassificationRequest = Field(description="Original classification request")
    tier: str = Field(description="Classification tier (e.g., 'Tier 1')")
    evidence: List[Dict[str, Any]] = Field(description="List of evidence dictionaries with AMP/ASCO/CAP codes")
    final_interpretation: str = Field(description="Final clinical interpretation text")
