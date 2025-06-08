from typing import List
from pydantic import BaseModel, Field
from core.api_models import Variant
from core.engine import ClassifiedVariant


class ClassificationBatchRequest(BaseModel):
    """Request model for batch variant classification."""
    variants: List[Variant] = Field(description="List of variants to classify")
    tumor_type: str = Field(description="OncoTree tumor type code")


class ClassificationResult(BaseModel):
    """Response model containing classification results."""
    classified_variants: List[ClassifiedVariant] = Field(description="List of classified variants with tiers")
    tumor_type: str = Field(description="OncoTree tumor type code used for classification")
    total_variants: int = Field(description="Total number of variants processed")
    pertinent_negatives: List[str] = Field(default=[], description="Pertinent negative findings")
