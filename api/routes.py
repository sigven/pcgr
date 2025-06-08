from typing import List
from fastapi import APIRouter, Depends, HTTPException
from sqlmodel import Session
from api.models import ClassificationBatchRequest, ClassificationResult
from api.database import get_session
from core.engine import assign_tier, ClassifiedVariant
from wrappers.vep_wrapper import get_vep_annotation

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
