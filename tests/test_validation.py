import pytest
from api.validators import VariantValidator
from core.api_models import ValidationSummary
from api.database import get_session


class TestVariantValidator:
    """Test comprehensive variant validation."""
    
    def test_snv_hgvs_validation(self):
        """Test SNV HGVS format validation."""
        session = next(get_session())
        try:
            validator = VariantValidator(session)
        
            result = validator.validate_batch(["NM_004333.4:c.1799T>A"], "MEL")
            assert result.summary["errors"] == 0
            
            result = validator.validate_batch(["NM_004333.4:c.1799>A"], "MEL")
            assert result.summary["errors"] == 1
            assert "Invalid HGVS" in result.errors[0].message
        finally:
            session.close()
    
    def test_genomic_coordinate_validation(self):
        """Test genomic coordinate validation."""
        session = next(get_session())
        try:
            validator = VariantValidator(session)
            
            result = validator.validate_batch(["chr7:140753336A>T"], "MEL")
            assert result.summary["errors"] == 0
            
            result = validator.validate_batch(["chr7:140753336A>"], "MEL")
            assert result.summary["errors"] == 1
        finally:
            session.close()
    
    def test_cnv_validation(self):
        """Test CNV validation."""
        session = next(get_session())
        try:
            validator = VariantValidator(session)
            
            result = validator.validate_batch(["chr17:4123-4127:copy_number=3.5"], "BRCA")
            assert result.summary["errors"] == 0
            
            result = validator.validate_batch(["chr17:4123-4127"], "BRCA")
            assert result.summary["errors"] == 1
            assert "Invalid CNV format" in result.errors[0].message
        finally:
            session.close()
    
    def test_fusion_validation(self):
        """Test fusion validation."""
        session = next(get_session())
        try:
            validator = VariantValidator(session)
            
            result = validator.validate_batch(["NM_000059.3:exon1-NM_007294.3:exon2"], "LAML")
            assert result.summary["errors"] == 0
            
            result = validator.validate_batch(["TMPRSS2-ERG"], "PRAD")
            assert result.summary["errors"] == 1
            assert "requires transcript IDs" in result.errors[0].message
        finally:
            session.close()
    
    def test_tumor_marker_validation(self):
        """Test tumor marker validation."""
        session = next(get_session())
        try:
            validator = VariantValidator(session)
            
            test_cases = [
                "MSI-H",
                "TMB:15.2",
                "HRD:45.0",
                "BRAF_expression:2.5"
            ]
            
            for case in test_cases:
                result = validator.validate_batch([case], "LUAD")
                assert result.summary["errors"] == 0, f"Failed for {case}"
        finally:
            session.close()
    
    def test_oncotree_validation(self):
        """Test OncoTree code validation."""
        session = next(get_session())
        try:
            validator = VariantValidator(session)
            
            result = validator.validate_batch(["chr7:140753336A>T"], "INVALID_CODE")
            assert result.summary["errors"] == 1
            assert "Invalid OncoTree code" in result.errors[0].message
        finally:
            session.close()
    
    def test_aggregate_error_reporting(self):
        """Test aggregate error reporting."""
        session = next(get_session())
        try:
            validator = VariantValidator(session)
            
            mixed_input = [
                "NM_004333.4:c.1799T>A",
                "invalid_format",
                "chr7:140753336A>T",
                "TMPRSS2-ERG"
            ]
            
            result = validator.validate_batch(mixed_input, "MEL")
            
            assert result.summary["total_variants"] == 4
            assert result.summary["parsed"] == 2
            assert result.summary["errors"] == 2
            assert len(result.errors) == 2
        finally:
            session.close()
