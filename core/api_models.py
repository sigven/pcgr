from typing import Union, List, Dict, Any, Optional
from pydantic import BaseModel, Field


class Variant(BaseModel):
    """Unified variant model supporting SNV, CNV, and SV types."""
    variant_type: str = Field(description="Type of variant: SNV, CNV, or SV")
    chromosome: str = Field(description="Chromosome (e.g., 'chr1', '1')")
    
    position: Optional[int] = Field(default=None, description="Genomic position for SNVs")
    reference: Optional[str] = Field(default=None, description="Reference allele for SNVs")
    alternate: Optional[str] = Field(default=None, description="Alternate allele for SNVs")
    
    start_position: Optional[int] = Field(default=None, description="Start position for CNVs/SVs")
    end_position: Optional[int] = Field(default=None, description="End position for CNVs")
    copy_number: Optional[int] = Field(default=None, description="Copy number for CNVs")
    
    sv_type: Optional[str] = Field(default=None, description="Structural variant type (DEL, DUP, INV, etc.)")
    
    reference_allele: Optional[str] = Field(default=None, description="Legacy reference allele field")
    alternate_allele: Optional[str] = Field(default=None, description="Legacy alternate allele field")
    gene: Optional[str] = Field(default=None, description="Gene symbol for CNVs")
    copy_number_change: Optional[str] = Field(default=None, description="Legacy CNV change description")
    genes: Optional[List[str]] = Field(default=None, description="List of genes for SVs")


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
