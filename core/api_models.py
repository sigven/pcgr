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
    """SNV/Indel with HGVS or genomic coordinate support."""
    transcript_id: Optional[str] = Field(default=None, description="RefSeq/Ensembl transcript ID (e.g., NM_000059.3)")
    hgvs_c: Optional[str] = Field(default=None, description="HGVS coding notation (e.g., c.725G>A)")
    hgvs_p: Optional[str] = Field(default=None, description="HGVS protein notation (e.g., p.Val600Glu)")
    chromosome: Optional[str] = Field(default=None, description="Chromosome (e.g., 'chr7')")
    position: Optional[int] = Field(default=None, description="Genomic position")
    reference: Optional[str] = Field(default=None, description="Reference allele")
    alternate: Optional[str] = Field(default=None, description="Alternate allele")
    build: str = Field(default="GRCh38", description="Genome build (GRCh37/GRCh38)")
    strand: Optional[str] = Field(default=None, description="Strand orientation")


class CNVInput(BaseModel):
    """Copy Number Variation."""
    chromosome: str = Field(description="Chromosome")
    start: int = Field(ge=1, description="Start position")
    end: int = Field(ge=1, description="End position")
    copy_number: float = Field(ge=0, description="Copy number value")
    build: str = Field(default="GRCh38", description="Genome build")


class SVInput(BaseModel):
    """Structural Variant."""
    breakends: List[Dict[str, Any]] = Field(description="List of breakend objects with chr, pos, strand, mate_id")
    svtype: str = Field(description="SV type: DEL, DUP, INV, TRA, COMPLEX")
    build: str = Field(default="GRCh38", description="Genome build")


class FusionInput(BaseModel):
    """Gene Fusion."""
    transcript5p: str = Field(description="5' transcript ID")
    junction5p_exon: str = Field(description="5' junction exon")
    transcript3p: str = Field(description="3' transcript ID")
    junction3p_exon: str = Field(description="3' junction exon")
    strand5p: Optional[str] = Field(default=None, description="5' strand")
    strand3p: Optional[str] = Field(default=None, description="3' strand")
    build: str = Field(default="GRCh38", description="Genome build")


class TumorMarkerInput(BaseModel):
    """Tumor markers and scores."""
    msi_status: Optional[str] = Field(default=None, description="MSI-H, MSI-L, MSS")
    msi_score: Optional[float] = Field(default=None, description="MSI fraction")
    tmb_score: Optional[float] = Field(default=None, description="TMB mutations/Mb")
    hrd_score: Optional[float] = Field(default=None, description="HRD score")
    expression_data: Optional[Dict[str, float]] = Field(default=None, description="Gene expression values")


class ValidationError(BaseModel):
    """Individual validation error."""
    line: int = Field(description="Line number in input")
    input: str = Field(description="Original input text")
    message: str = Field(description="Error description")


class ValidationSummary(BaseModel):
    """Aggregate validation report."""
    summary: Dict[str, int] = Field(description="Total, parsed, errors counts")
    errors: List[ValidationError] = Field(description="List of validation errors")


class QCFlags(BaseModel):
    """Quality control flags for variant analysis."""
    qns: Optional[bool] = Field(default=None, description="Quantity Not Sufficient")
    low_read_depth: Optional[bool] = Field(default=None, description="Low read depth detected")
    failed_to_amplify: Optional[bool] = Field(default=None, description="Failed to amplify")
    poor_quality: Optional[bool] = Field(default=None, description="Poor sample quality")
    contamination: Optional[bool] = Field(default=None, description="Sample contamination detected")
    notes: Optional[str] = Field(default=None, description="Additional QC notes")


class AttendingInterpretation(BaseModel):
    """Attending physician interpretation."""
    attending_name: str = Field(description="Name of attending physician")
    interpretation: Optional[str] = Field(default=None, description="Attending's clinical interpretation")
    submitted_before_tiering: Optional[bool] = Field(default=None, description="Whether interpretation was submitted before seeing tier assignments")
    confidence_level: Optional[str] = Field(default=None, description="Confidence level: high, medium, low")
    clinical_significance: Optional[str] = Field(default=None, description="Clinical significance assessment")


class ClassificationRequest(BaseModel):
    """Main model for an API classification request."""
    variant: Variant = Field(description="Variant to classify (SNV, CNV, or SV)")
    oncotree_code: str = Field(description="OncoTree disease code")
    attending_interpretation: Optional[AttendingInterpretation] = Field(default=None, description="Optional attending physician interpretation")
    qc_flags: Optional[QCFlags] = Field(default=None, description="Optional quality control flags")


class ClassificationReport(BaseModel):
    """Final report object returned by the API."""
    request: ClassificationRequest = Field(description="Original classification request")
    tier: str = Field(description="Classification tier (e.g., 'Tier 1')")
    evidence: List[Dict[str, Any]] = Field(description="List of evidence dictionaries with AMP/ASCO/CAP codes")
    final_interpretation: str = Field(description="Final clinical interpretation text")
    attending_review: Optional[Dict[str, Any]] = Field(default=None, description="Attending physician review and interpretation")
    qc_summary: Optional[Dict[str, Any]] = Field(default=None, description="Quality control summary and flags")
