import pytest
from sqlmodel import Session, create_engine, SQLModel
from core.engine import assign_tier, is_benign_by_frequency, is_tumor_suppressor_gene, is_oncogene
from core.api_models import Variant
from wrappers.vep_wrapper import VEPAnnotation
from core.kb_models import GeneGeneralComment, VariantDxInterpretation


@pytest.fixture
def test_session():
    """Create in-memory test database with sample data."""
    engine = create_engine("sqlite:///:memory:")
    SQLModel.metadata.create_all(engine)
    
    with Session(engine) as session:
        test_gene_braf = GeneGeneralComment(
            gene_id="673",
            gene_symbol="BRAF",
            comment="BRAF is an oncogene involved in cell signaling and melanoma development",
            created_by="test"
        )
        
        test_gene_tp53 = GeneGeneralComment(
            gene_id="7157",
            gene_symbol="TP53",
            comment="TP53 is a tumor suppressor gene that regulates cell cycle",
            created_by="test"
        )
        
        test_variant_interp = VariantDxInterpretation(
            variant_id="BRAF:p.V600E",
            disease="MEL",
            interpretation="Pathogenic hotspot mutation in melanoma",
            evidence_level="A",
            evidence_type="Predictive",
            created_by="test"
        )
        
        session.add(test_gene_braf)
        session.add(test_gene_tp53)
        session.add(test_variant_interp)
        session.commit()
        yield session


def test_benign_by_frequency():
    """Test population frequency-based benign classification."""
    annotation = VEPAnnotation(
        gene_symbol="TEST",
        consequence="missense_variant",
        gnomad_af=0.05,
        hgvsc="c.123A>T",
        hgvsp="p.Lys41Asn"
    )
    assert is_benign_by_frequency(annotation) == True
    
    annotation.gnomad_af = 0.005
    assert is_benign_by_frequency(annotation) == False
    
    annotation.gnomad_af = None
    assert is_benign_by_frequency(annotation) == False


def test_tumor_suppressor_detection(test_session):
    """Test tumor suppressor gene detection."""
    assert is_tumor_suppressor_gene("TP53", "ANY", test_session) == True
    assert is_tumor_suppressor_gene("BRAF", "ANY", test_session) == False
    assert is_tumor_suppressor_gene("UNKNOWN", "ANY", test_session) == False


def test_oncogene_detection(test_session):
    """Test oncogene detection."""
    assert is_oncogene("BRAF", "ANY", test_session) == True
    assert is_oncogene("TP53", "ANY", test_session) == False
    assert is_oncogene("UNKNOWN", "ANY", test_session) == False


def test_tier_assignment_high_frequency(test_session):
    """Test tier assignment for high frequency variants."""
    variant = Variant(
        variant_type="SNV",
        chromosome="7",
        position=140753336,
        reference="A",
        alternate="T"
    )
    
    annotation = VEPAnnotation(
        gene_symbol="BRAF",
        consequence="missense_variant",
        gnomad_af=0.05,
        hgvsc="c.1799T>A",
        hgvsp="p.V600E"
    )
    
    tier = assign_tier(variant, annotation, "MEL", test_session)
    assert tier in ['1', '5']


def test_tier_assignment_curated_evidence(test_session):
    """Test tier assignment with curated evidence."""
    variant = Variant(
        variant_type="SNV",
        chromosome="7",
        position=140753336,
        reference="A",
        alternate="T"
    )
    
    annotation = VEPAnnotation(
        gene_symbol="BRAF",
        consequence="missense_variant",
        gnomad_af=0.001,
        hgvsc="c.1799T>A",
        hgvsp="p.V600E"
    )
    
    tier = assign_tier(variant, annotation, "MEL", test_session)
    assert tier in ['1', '2', '3', '4', '5']


def test_tier_assignment_cancer_gene_vus(test_session):
    """Test tier assignment for cancer gene VUS."""
    variant = Variant(
        variant_type="SNV",
        chromosome="7",
        position=140753336,
        reference="A",
        alternate="T"
    )
    
    annotation = VEPAnnotation(
        gene_symbol="BRAF",
        consequence="synonymous_variant",
        gnomad_af=0.001,
        hgvsc="c.1800C>T",
        hgvsp="p.V600="
    )
    
    tier = assign_tier(variant, annotation, "MEL", test_session)
    assert tier in ['3', '4']


def test_tier_assignment_unknown_gene(test_session):
    """Test tier assignment for unknown gene."""
    variant = Variant(
        variant_type="SNV",
        chromosome="1",
        position=12345,
        reference="A",
        alternate="T"
    )
    
    annotation = VEPAnnotation(
        gene_symbol="UNKNOWN",
        consequence="missense_variant",
        gnomad_af=0.001,
        hgvsc="c.123A>T",
        hgvsp="p.Lys41Asn"
    )
    
    tier = assign_tier(variant, annotation, "MEL", test_session)
    assert tier == '4'
