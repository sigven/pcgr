import pytest
from fastapi.testclient import TestClient
from sqlmodel import Session, create_engine, SQLModel
from sqlalchemy.pool import StaticPool
from api.main import app
from api.database import get_session
from core.kb_models import GeneGeneralComment, VariantDxInterpretation


@pytest.fixture
def test_client():
    """Create test client with in-memory database."""
    engine = create_engine(
        "sqlite:///:memory:",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool
    )
    SQLModel.metadata.create_all(engine)
    
    def get_test_session():
        with Session(engine) as session:
            yield session
    
    app.dependency_overrides[get_session] = get_test_session
    
    with Session(engine) as session:
        test_gene = GeneGeneralComment(
            gene_id="673",
            gene_symbol="BRAF",
            comment="BRAF is an oncogene",
            created_by="test"
        )
        
        test_variant = VariantDxInterpretation(
            variant_id="BRAF:p.V600E",
            disease="MEL",
            interpretation="Pathogenic mutation",
            evidence_level="A",
            created_by="test"
        )
        
        session.add(test_gene)
        session.add(test_variant)
        session.commit()
    
    client = TestClient(app)
    yield client
    
    app.dependency_overrides.clear()


def test_classify_endpoint(test_client):
    """Test main classification endpoint."""
    request_data = {
        "variants": [{
            "variant_type": "SNV",
            "chromosome": "7",
            "position": 140753336,
            "reference": "A",
            "alternate": "T"
        }],
        "tumor_type": "MEL"
    }
    
    response = test_client.post("/v1/classify", json=request_data)
    assert response.status_code == 200
    data = response.json()
    assert "classified_variants" in data
    assert len(data["classified_variants"]) == 1
    assert data["tumor_type"] == "MEL"


def test_validate_variant_endpoint_valid_snv(test_client):
    """Test variant validation endpoint with valid SNV."""
    valid_variant = {
        "variant_type": "SNV",
        "chromosome": "7",
        "position": 140753336,
        "reference": "A",
        "alternate": "T"
    }
    
    response = test_client.post("/v1/validate_variant", json=valid_variant)
    assert response.status_code == 200
    assert response.json()["valid"] == True


def test_validate_variant_endpoint_valid_cnv(test_client):
    """Test variant validation endpoint with valid CNV."""
    valid_variant = {
        "variant_type": "CNV",
        "chromosome": "7",
        "start_position": 140753336,
        "end_position": 140753400,
        "copy_number": 4
    }
    
    response = test_client.post("/v1/validate_variant", json=valid_variant)
    assert response.status_code == 200
    assert response.json()["valid"] == True


def test_validate_variant_endpoint_valid_sv(test_client):
    """Test variant validation endpoint with valid SV."""
    valid_variant = {
        "variant_type": "SV",
        "chromosome": "7",
        "start_position": 140753336,
        "sv_type": "DEL"
    }
    
    response = test_client.post("/v1/validate_variant", json=valid_variant)
    assert response.status_code == 200
    assert response.json()["valid"] == True


def test_validate_variant_endpoint_invalid_type(test_client):
    """Test variant validation endpoint with invalid type."""
    invalid_variant = {
        "variant_type": "INVALID",
        "chromosome": "7",
        "position": 140753336
    }
    
    response = test_client.post("/v1/validate_variant", json=invalid_variant)
    assert response.status_code == 200
    assert response.json()["valid"] == False
    assert "Invalid variant type" in response.json()["message"]


def test_validate_variant_endpoint_missing_snv_fields(test_client):
    """Test variant validation endpoint with missing SNV fields."""
    invalid_variant = {
        "variant_type": "SNV",
        "chromosome": "7"
    }
    
    response = test_client.post("/v1/validate_variant", json=invalid_variant)
    assert response.status_code == 200
    assert response.json()["valid"] == False
    assert "missing required fields" in response.json()["message"]


def test_gene_interpretations_endpoint(test_client):
    """Test gene interpretations lookup."""
    response = test_client.get("/v1/interpretations/gene/BRAF")
    assert response.status_code == 200
    data = response.json()
    assert data["gene_symbol"] == "BRAF"
    assert "general_comments" in data
    assert "variant_interpretations" in data
    assert "pertinent_negatives" in data
    assert len(data["general_comments"]) > 0
    assert len(data["variant_interpretations"]) > 0


def test_gene_interpretations_endpoint_unknown_gene(test_client):
    """Test gene interpretations lookup for unknown gene."""
    response = test_client.get("/v1/interpretations/gene/UNKNOWN")
    assert response.status_code == 200
    data = response.json()
    assert data["gene_symbol"] == "UNKNOWN"
    assert len(data["general_comments"]) == 0
    assert len(data["variant_interpretations"]) == 0


def test_report_endpoint(test_client):
    """Test report retrieval endpoint."""
    response = test_client.get("/v1/report/test-analysis-123")
    assert response.status_code == 200
    data = response.json()
    assert data["analysis_id"] == "test-analysis-123"
    assert data["status"] == "completed"
    assert "report_url" in data
    assert "generated_at" in data


def test_classify_endpoint_multiple_variants(test_client):
    """Test classification with multiple variants."""
    request_data = {
        "variants": [
            {
                "variant_type": "SNV",
                "chromosome": "7",
                "position": 140753336,
                "reference": "A",
                "alternate": "T"
            },
            {
                "variant_type": "CNV",
                "chromosome": "17",
                "start_position": 7565097,
                "end_position": 7590856,
                "copy_number": 0
            }
        ],
        "tumor_type": "MEL"
    }
    
    response = test_client.post("/v1/classify", json=request_data)
    assert response.status_code == 200
    data = response.json()
    assert len(data["classified_variants"]) == 2
    assert data["total_variants"] == 2
