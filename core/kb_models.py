from datetime import datetime
from typing import List, Optional
from sqlmodel import SQLModel, Field


class GeneGeneralComment(SQLModel, table=True):
    """General, non-disease-specific information about a gene.
    
    Based on: resources/gene_text/general_gene_comments.tsv
    Schema: GeneID, Symbol, Summary
    """
    id: Optional[int] = Field(default=None, primary_key=True)
    gene_id: str = Field(description="Gene identifier (e.g., Entrez Gene ID)")
    gene_symbol: str = Field(description="Gene symbol (e.g., 'BRAF')")
    comment: str = Field(description="General gene summary/description")
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    created_by: str


class VariantGeneralComment(SQLModel, table=True):
    """General, non-disease-specific information about a specific variant.
    
    Based on: resources/variant_interps/general_variant_comments.tsv
    Schema: variant_id, gene, variant, blurb, chromosome, start, stop, variant_types, hgvs_descriptions, etc.
    """
    id: Optional[int] = Field(default=None, primary_key=True)
    variant_id: str = Field(description="Standardized variant identifier")
    gene: str = Field(description="Gene symbol")
    variant: str = Field(description="Variant description")
    comment: str = Field(description="General variant commentary", alias="blurb")
    chromosome: Optional[str] = Field(default=None, description="Chromosome")
    start: Optional[int] = Field(default=None, description="Start position")
    stop: Optional[int] = Field(default=None, description="Stop position")
    variant_types: Optional[str] = Field(default=None, description="Variant type classification")
    hgvs_descriptions: Optional[str] = Field(default=None, description="HGVS notation")
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    created_by: str


class OncoTreeMapping(SQLModel, table=True):
    """Maps disease names to OncoTree codes.
    
    Note: This structure is inferred from the variant_dx_interpretations.tsv 
    which contains disease and doid columns that would map to OncoTree.
    """
    id: Optional[int] = Field(default=None, primary_key=True)
    disease_name: str = Field(description="Human-readable disease name (e.g., 'Melanoma')")
    oncotree_code: str = Field(description="OncoTree code (e.g., 'MEL')")
    doid: Optional[str] = Field(default=None, description="Disease Ontology ID")
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    created_by: str


class VariantDxInterpretation(SQLModel, table=True):
    """Core model for disease-specific clinical interpretations of a variant.
    
    Based on: resources/variant_interps/variant_dx_interpretations.tsv
    Schema: molecular_profile, molecular_profile_id, disease, doid, phenotypes, therapies, 
    therapy_interaction_type, evidence_type, evidence_direction, evidence_level, significance, 
    evidence_statement, citation_id, source_type, asco_abstract_id, citation, nct_ids, rating, 
    evidence_status, evidence_id, variant_origin, last_review_date, evidence_civic_url, 
    molecular_profile_civic_url, is_flagged
    """
    id: Optional[int] = Field(default=None, primary_key=True)
    variant_id: str = Field(description="Variant identifier", alias="molecular_profile")
    molecular_profile_id: Optional[str] = Field(default=None, description="CIViC molecular profile ID")
    disease: str = Field(description="Disease name (OncoTree code when available)")
    doid: Optional[str] = Field(default=None, description="Disease Ontology ID")
    interpretation: str = Field(description="Clinical interpretation text", alias="evidence_statement")
    evidence_type: Optional[str] = Field(default=None, description="Type of evidence (Diagnostic, Predictive, etc.)")
    evidence_direction: Optional[str] = Field(default=None, description="Evidence direction (Supports, Does Not Support)")
    evidence_level: Optional[str] = Field(default=None, description="Evidence level (A, B, C, D, E)")
    significance: Optional[str] = Field(default=None, description="Clinical significance")
    references: str = Field(default="", description="Citation references (JSON array)")
    citation_id: Optional[str] = Field(default=None, description="Primary citation ID")
    citation: Optional[str] = Field(default=None, description="Citation text")
    nct_ids: Optional[str] = Field(default=None, description="Clinical trial NCT IDs")
    evidence_status: Optional[str] = Field(default=None, description="Evidence status (accepted, rejected, etc.)")
    last_review_date: Optional[datetime] = Field(default=None, description="Last review date")
    civic_url: Optional[str] = Field(default=None, description="CIViC evidence URL", alias="evidence_civic_url")
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    created_by: str


class DrugAssociation(SQLModel, table=True):
    """Links variants and diseases to therapies.
    
    Based on: resources/variant_interps/variant_dx_interpretations.tsv (therapies column)
    and resources/Biomarkers/Biomarkers.tsv for additional drug associations
    """
    id: Optional[int] = Field(default=None, primary_key=True)
    variant_id: str = Field(description="Variant identifier")
    disease: str = Field(description="Disease (OncoTree code)")
    drug_name: str = Field(description="Drug/therapy name")
    response_type: str = Field(description="Response type (Sensitivity, Resistance, etc.)")
    references: str = Field(default="", description="Supporting citations (JSON array)")
    evidence_level: Optional[str] = Field(default=None, description="Evidence level")
    therapy_interaction_type: Optional[str] = Field(default=None, description="Interaction type")
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    created_by: str


class PertinentNegative(SQLModel, table=True):
    """Genes important to report when no pathogenic variants are found in specific disease context.
    
    Based on: resources/pertinent_negatives/Pertinent_Negatives.tsv
    Schema: Type, Tissue, Pertinent Negatives, Notes
    """
    id: Optional[int] = Field(default=None, primary_key=True)
    gene_symbol: str = Field(description="Gene symbol to report as negative")
    disease: str = Field(description="Disease context (tissue type or OncoTree code)", alias="Tissue")
    comment: str = Field(description="Explanatory comment", alias="Notes")
    pertinent_negatives_text: Optional[str] = Field(default=None, description="Raw pertinent negatives text", alias="Pertinent Negatives")
    negative_type: Optional[str] = Field(default=None, description="Type of pertinent negative", alias="Type")
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    created_by: str
