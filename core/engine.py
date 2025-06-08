from typing import Optional
from sqlmodel import Session, select
from pydantic import BaseModel, Field

from core.api_models import Variant
from core.kb_models import VariantDxInterpretation, GeneGeneralComment
from wrappers.vep_wrapper import VEPAnnotation


class ClassifiedVariant(BaseModel):
    """Model for a variant with its annotation and assigned tier."""
    variant: Variant = Field(description="Original variant input")
    annotation: VEPAnnotation = Field(description="VEP annotation data")
    tier: str = Field(description="Assigned classification tier (1-5)")
    evidence_summary: str = Field(description="Summary of evidence used for classification")


def assign_tier(variant: Variant, annotation: VEPAnnotation, tumor_type: str, db_session: Session) -> str:
    """
    Core tiering engine that classifies variants using prioritized evidence-based rules.
    
    Args:
        variant: The input variant (SNV, CNV, or SV)
        annotation: VEP annotation data including consequences and population frequencies
        tumor_type: OncoTree tumor type code
        db_session: SQLModel database session for knowledge base queries
        
    Returns:
        Tier string ('1', '2', '3', '4', '5')
    """
    
    tier_from_kb = check_curated_evidence(annotation, tumor_type, db_session)
    if tier_from_kb:
        return tier_from_kb
    
    if is_benign_by_frequency(annotation):
        return '5'
    
    tier_from_oncogenic = check_oncogenic_evidence(annotation, tumor_type, db_session)
    if tier_from_oncogenic:
        return tier_from_oncogenic
    
    if is_cancer_gene_vus(annotation, db_session):
        return '3'
    
    return '4'


def check_curated_evidence(annotation: VEPAnnotation, tumor_type: str, db_session: Session) -> Optional[str]:
    """
    Priority 1: Query VariantDxInterpretation table for exact matches.
    
    Returns tier from database record if exact match found for hgvsc/hgvsp and tumor type.
    """
    if not annotation.hgvsc and not annotation.hgvsp:
        return None
    
    query = select(VariantDxInterpretation).where(
        VariantDxInterpretation.disease == tumor_type
    )
    
    if annotation.hgvsc:
        hgvsc_query = query.where(VariantDxInterpretation.variant_id.contains(annotation.hgvsc))
        result = db_session.exec(hgvsc_query).first()
        if result:
            return extract_tier_from_evidence_level(result.evidence_level)
    
    if annotation.hgvsp:
        hgvsp_query = query.where(VariantDxInterpretation.variant_id.contains(annotation.hgvsp))
        result = db_session.exec(hgvsp_query).first()
        if result:
            return extract_tier_from_evidence_level(result.evidence_level)
    
    return None


def is_benign_by_frequency(annotation: VEPAnnotation) -> bool:
    """
    Priority 2: Check if variant is benign based on population frequency.
    
    Returns True if gnomAD AF > 1% (somatic variant threshold).
    """
    if annotation.gnomad_af is None:
        return False
    
    return annotation.gnomad_af > 0.01


def check_oncogenic_evidence(annotation: VEPAnnotation, tumor_type: str, db_session: Session) -> Optional[str]:
    """
    Priority 3: Check for strong oncogenic evidence using rule-based logic.
    
    Returns '1' or '2' for strong evidence, None otherwise.
    """
    lof_consequences = {
        'frameshift_variant',
        'stop_gained', 
        'stop_lost',
        'start_lost',
        'splice_donor_variant',
        'splice_acceptor_variant'
    }
    
    if annotation.consequence in lof_consequences:
        if is_tumor_suppressor_gene(annotation.gene_symbol, tumor_type, db_session):
            return '1'  # Strong evidence for pathogenicity
    
    if is_oncogenic_hotspot(annotation, db_session):
        return '1'  # Strong evidence for pathogenicity
    
    if annotation.consequence == 'missense_variant':
        if is_oncogene(annotation.gene_symbol, tumor_type, db_session):
            return '2'  # Moderate evidence for pathogenicity
    
    return None


def is_cancer_gene_vus(annotation: VEPAnnotation, db_session: Session) -> bool:
    """
    Priority 4: Check if variant is in a known cancer gene but doesn't meet strong criteria.
    
    Returns True if gene exists in knowledge base (indicating cancer relevance).
    """
    query = select(GeneGeneralComment).where(
        GeneGeneralComment.gene_symbol == annotation.gene_symbol
    )
    result = db_session.exec(query).first()
    return result is not None


def is_tumor_suppressor_gene(gene_symbol: str, tumor_type: str, db_session: Session) -> bool:
    """Check if gene is a known tumor suppressor for the given tumor type."""
    query = select(GeneGeneralComment).where(
        GeneGeneralComment.gene_symbol == gene_symbol
    )
    result = db_session.exec(query).first()
    
    if result and result.comment:
        tsg_keywords = ['tumor suppressor', 'tumour suppressor', 'TSG', 'loss of function']
        comment_lower = result.comment.lower()
        return any(keyword in comment_lower for keyword in tsg_keywords)
    
    return False


def is_oncogene(gene_symbol: str, tumor_type: str, db_session: Session) -> bool:
    """Check if gene is a known oncogene for the given tumor type."""
    query = select(GeneGeneralComment).where(
        GeneGeneralComment.gene_symbol == gene_symbol
    )
    result = db_session.exec(query).first()
    
    if result and result.comment:
        oncogene_keywords = ['oncogene', 'gain of function', 'activating', 'driver']
        comment_lower = result.comment.lower()
        return any(keyword in comment_lower for keyword in oncogene_keywords)
    
    return False


def is_oncogenic_hotspot(annotation: VEPAnnotation, db_session: Session) -> bool:
    """Check if variant matches a known oncogenic hotspot mutation."""
    if not annotation.hgvsp:
        return False
    
    query = select(VariantDxInterpretation).where(
        VariantDxInterpretation.variant_id.contains(annotation.hgvsp)
    )
    results = db_session.exec(query).all()
    
    for result in results:
        if result.interpretation and result.evidence_level:
            interp_lower = result.interpretation.lower()
            hotspot_keywords = ['hotspot', 'pathogenic', 'oncogenic', 'activating']
            if any(keyword in interp_lower for keyword in hotspot_keywords):
                return True
    
    return False


def extract_tier_from_evidence_level(evidence_level: Optional[str]) -> str:
    """Convert evidence level to tier classification."""
    if not evidence_level:
        return '3'  # Default to VUS if no evidence level
    
    level = evidence_level.upper()
    
    if level in ['A', 'LEVEL_A', '1']:
        return '1'  # Strong clinical significance
    elif level in ['B', 'LEVEL_B', '2']:
        return '2'  # Moderate clinical significance  
    elif level in ['C', 'LEVEL_C', '3']:
        return '3'  # Uncertain significance
    elif level in ['D', 'LEVEL_D', '4']:
        return '4'  # Likely benign
    elif level in ['E', 'LEVEL_E', '5']:
        return '5'  # Benign
    else:
        return '3'  # Default to VUS for unknown levels
