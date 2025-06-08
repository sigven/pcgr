from typing import List
from fastapi import APIRouter, Depends, HTTPException
from sqlmodel import Session, select
from api.models import ClassificationBatchRequest, ClassificationResult
from api.database import get_session
from core.engine import assign_tier, ClassifiedVariant
from core.api_models import Variant
from core.kb_models import GeneGeneralComment, VariantDxInterpretation, PertinentNegative
from wrappers.vep_wrapper import get_vep_annotation
from datetime import datetime

router = APIRouter()


@router.post("/classify", response_model=ClassificationResult)
async def classify_variants(
    request: ClassificationBatchRequest,
    db_session: Session = Depends(get_session)
) -> ClassificationResult:
    """
    Main endpoint for variant classification.
    
    Accepts a list of variants and returns complete classification report
    with tier assignments and evidence summaries.
    """
    classified_variants = []
    
    for variant in request.variants:
        try:
            annotation = get_vep_annotation(variant)
            
            tier = assign_tier(
                variant=variant,
                annotation=annotation,
                tumor_type=request.tumor_type,
                db_session=db_session
            )
            
            evidence_summary = generate_evidence_summary(tier, annotation)
            
            classified_variant = ClassifiedVariant(
                variant=variant,
                annotation=annotation,
                tier=tier,
                evidence_summary=evidence_summary
            )
            
            classified_variants.append(classified_variant)
            
        except Exception as e:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to classify variant: {str(e)}"
            )
    
    return ClassificationResult(
        classified_variants=classified_variants,
        tumor_type=request.tumor_type,
        total_variants=len(classified_variants),
        pertinent_negatives=[]
    )


@router.post("/validate_variant")
async def validate_variant(variant: Variant) -> dict:
    """Validate variant input format before classification."""
    try:
        if variant.variant_type not in ['SNV', 'CNV', 'SV']:
            raise ValueError("Invalid variant type")
        
        if variant.variant_type == 'SNV':
            if not all([variant.chromosome, variant.position, variant.reference, variant.alternate]):
                raise ValueError("SNV missing required fields: chromosome, position, reference, alternate")
        
        if variant.variant_type == 'CNV':
            if not all([variant.chromosome, variant.start_position, variant.end_position, variant.copy_number]):
                raise ValueError("CNV missing required fields: chromosome, start_position, end_position, copy_number")
        
        if variant.variant_type == 'SV':
            if not all([variant.chromosome, variant.start_position, variant.sv_type]):
                raise ValueError("SV missing required fields: chromosome, start_position, sv_type")
        
        return {"valid": True, "message": "Variant format is valid"}
    except Exception as e:
        return {"valid": False, "message": str(e)}


@router.get("/interpretations/gene/{gene_symbol}")
async def get_gene_interpretations(
    gene_symbol: str,
    db_session: Session = Depends(get_session)
) -> dict:
    """Retrieve all interpretations for a specific gene."""
    gene_comments = db_session.exec(
        select(GeneGeneralComment).where(GeneGeneralComment.gene_symbol == gene_symbol)
    ).all()
    
    variant_interps = db_session.exec(
        select(VariantDxInterpretation).where(VariantDxInterpretation.variant_id.contains(gene_symbol))
    ).all()
    
    pertinent_negatives = db_session.exec(
        select(PertinentNegative).where(PertinentNegative.gene_symbol == gene_symbol)
    ).all()
    
    return {
        "gene_symbol": gene_symbol,
        "general_comments": [{"comment": c.comment, "gene_id": c.gene_id} for c in gene_comments],
        "variant_interpretations": [
            {
                "disease": v.disease,
                "interpretation": v.interpretation,
                "evidence_level": v.evidence_level,
                "evidence_type": v.evidence_type
            } for v in variant_interps
        ],
        "pertinent_negatives": [
            {
                "disease": p.disease,
                "comment": p.comment
            } for p in pertinent_negatives
        ]
    }


@router.get("/report/{analysis_id}")
async def get_report(analysis_id: str) -> dict:
    """Retrieve generated report by analysis ID."""
    return {
        "analysis_id": analysis_id,
        "status": "completed",
        "report_url": f"/reports/{analysis_id}.html",
        "generated_at": datetime.utcnow().isoformat(),
        "message": "Report generation system not yet implemented"
    }


def generate_evidence_summary(tier: str, annotation) -> str:
    """Generate human-readable evidence summary for classification."""
    tier_descriptions = {
        '1': f"Strong clinical significance - Gene: {annotation.gene_symbol}, Consequence: {annotation.consequence}",
        '2': f"Moderate clinical significance - Gene: {annotation.gene_symbol}, Consequence: {annotation.consequence}",
        '3': f"Variant of uncertain significance - Gene: {annotation.gene_symbol}, Consequence: {annotation.consequence}",
        '4': f"Likely benign - Gene: {annotation.gene_symbol}, Consequence: {annotation.consequence}",
        '5': f"Benign (population frequency: {annotation.gnomad_af}) - Gene: {annotation.gene_symbol}"
    }
    
    return tier_descriptions.get(tier, f"Classification tier {tier} - Gene: {annotation.gene_symbol}")
