from typing import List, Dict, Any, Tuple, Optional, Union
from pydantic import ValidationError as PydanticValidationError
import re
import os
import pandas as pd
from core.api_models import SNVInput, CNVInput, SVInput, FusionInput, TumorMarkerInput, ValidationError, ValidationSummary, QCFlags, AttendingInterpretation
from core.kb_models import OncoTreeMapping
from sqlmodel import Session, select


class VariantValidator:
    """Comprehensive variant validation with transcript lookup."""
    
    def __init__(self, db_session: Session):
        self.db_session = db_session
        self.transcript_lookup = self._build_transcript_lookup()
        self.oncotree_codes = self._load_oncotree_codes()
        self.qc_flags_config = self._load_qc_flags()
    
    def _build_transcript_lookup(self) -> Dict[str, Dict[str, str]]:
        """Build transcript ID lookup table from knowledge base."""
        transcript_map = {}
        
        common_transcripts = {
            "NM_000059.3": {"gene": "BRCA2", "type": "RefSeq"},
            "NM_007294.3": {"gene": "BRCA1", "type": "RefSeq"},
            "NM_004333.4": {"gene": "BRAF", "type": "RefSeq"},
            "NM_000546.5": {"gene": "TP53", "type": "RefSeq"},
            "ENST00000288602": {"gene": "BRCA1", "type": "Ensembl"},
            "ENST00000544455": {"gene": "BRCA2", "type": "Ensembl"},
            "ENST00000288602.11": {"gene": "BRCA1", "type": "Ensembl"},
        }
        
        transcript_map.update(common_transcripts)
        return transcript_map
    
    def _load_oncotree_codes(self) -> set:
        """Load valid OncoTree codes from database."""
        try:
            oncotree_mappings = self.db_session.exec(select(OncoTreeMapping)).all()
            codes = {mapping.oncotree_code for mapping in oncotree_mappings}
            if not codes:
                default_codes = {
                    "MEL", "LUAD", "BRCA", "CRC", "PRAD", "BLCA", "HNSC", "KIRC", 
                    "LIHC", "THCA", "UCEC", "KIRP", "SARC", "LAML", "PAAD", "GBM",
                    "OV", "SKCM", "CESC", "TGCT", "UCS", "CHOL", "THYM", "ACC",
                    "MESO", "UVM", "DLBC", "KICH", "READ", "LGG", "PCPG", "STAD"
                }
                return default_codes
            return codes
        except Exception:
            default_codes = {
                "MEL", "LUAD", "BRCA", "CRC", "PRAD", "BLCA", "HNSC", "KIRC", 
                "LIHC", "THCA", "UCEC", "KIRP", "SARC", "LAML", "PAAD", "GBM",
                "OV", "SKCM", "CESC", "TGCT", "UCS", "CHOL", "THYM", "ACC",
                "MESO", "UVM", "DLBC", "KICH", "READ", "LGG", "PCPG", "STAD"
            }
            return default_codes
    
    def validate_batch(self, variant_lines: List[str], oncotree_code: str, attending_interpretation: Optional[AttendingInterpretation] = None, qc_flags: Optional[QCFlags] = None) -> ValidationSummary:
        """Validate batch of variant inputs with aggregate error reporting."""
        errors = []
        parsed_variants = []
        
        if oncotree_code not in self.oncotree_codes:
            errors.append(ValidationError(
                line=0,
                input=oncotree_code,
                message=f"Invalid OncoTree code. Supported codes: {sorted(list(self.oncotree_codes))[:10]}..."
            ))
        
        if attending_interpretation:
            if not attending_interpretation.attending_name or len(attending_interpretation.attending_name.strip()) == 0:
                errors.append(ValidationError(
                    line=0,
                    input="attending_interpretation",
                    message="Attending name is required when providing attending interpretation"
                ))
        
        if qc_flags:
            qc_errors = self._validate_qc_flags(qc_flags)
            errors.extend(qc_errors)
        
        for line_num, line in enumerate(variant_lines, 1):
            try:
                variant_type, parsed_variant = self._parse_variant_line(line)
                if parsed_variant:
                    parsed_variants.append((variant_type, parsed_variant))
            except Exception as e:
                errors.append(ValidationError(
                    line=line_num,
                    input=line,
                    message=str(e)
                ))
        
        return ValidationSummary(
            summary={
                "total_variants": len(variant_lines),
                "parsed": len(parsed_variants),
                "errors": len(errors)
            },
            errors=errors
        )
    
    def _parse_variant_line(self, line: str) -> Tuple[str, Optional[Any]]:
        """Parse single variant line and return type and validated object."""
        line = line.strip()
        
        if self._is_hgvs_format(line):
            return "SNV", self._validate_snv_hgvs(line)
        elif self._is_genomic_coordinate(line):
            return "SNV", self._validate_snv_genomic(line)
        elif self._is_cnv_format(line):
            return "CNV", self._validate_cnv(line)
        elif self._is_fusion_format(line):
            return "Fusion", self._validate_fusion(line)
        elif self._is_sv_format(line):
            return "SV", self._validate_sv(line)
        elif self._is_tumor_marker(line):
            return "TumorMarker", self._validate_tumor_marker(line)
        elif self._is_simple_fusion_format(line):
            raise ValueError("Fusion requires transcript IDs: e.g. 'NM_xxx:exon1-NM_yyy:exon2'")
        else:
            raise ValueError(f"Unrecognized variant format. Expected formats: transcript:hgvs, chr:pos:ref:alt, gene amplification, gene1-gene2 fusion, MSI-H, TMB:10.5")
    
    def _is_hgvs_format(self, line: str) -> bool:
        """Check if line matches HGVS format (transcript:c.notation or transcript:p.notation)."""
        return bool(re.match(r'^(NM_\d+\.\d+|ENST\d+):(c\.|p\.)', line))
    
    def _is_genomic_coordinate(self, line: str) -> bool:
        """Check if line matches genomic coordinate format."""
        return bool(re.match(r'^chr\w+:\d+[ATCG]>[ATCG]$', line, re.IGNORECASE))
    
    def _is_cnv_format(self, line: str) -> bool:
        """Check if line matches CNV format."""
        return bool(re.match(r'^chr\w+:\d+-\d+(:copy_number=[\d.]+)?$', line, re.IGNORECASE))
    
    def _is_fusion_format(self, line: str) -> bool:
        """Check if line matches fusion format."""
        return bool(re.match(r'^(NM_\d+\.\d+|ENST\d+):exon\w+-(NM_\d+\.\d+|ENST\d+):exon\w+$', line))
    
    def _is_simple_fusion_format(self, line: str) -> bool:
        """Check if line matches simple gene-gene fusion format (e.g., TMPRSS2-ERG)."""
        return bool(re.match(r'^[A-Z0-9]+[-_][A-Z0-9]+$', line))
    
    def _is_sv_format(self, line: str) -> bool:
        """Check if line matches SV format."""
        return bool(re.match(r'^(DEL|DUP|INV|TRA|COMPLEX):', line, re.IGNORECASE))
    
    def _is_tumor_marker(self, line: str) -> bool:
        """Check if line matches tumor marker format."""
        patterns = [
            r'^MSI-(H|L)$',
            r'^MSS$',
            r'^TMB:[\d.]+$',
            r'^HRD:[\d.]+$',
            r'^\w+_expression:[\d.]+$'
        ]
        return any(re.match(pattern, line, re.IGNORECASE) for pattern in patterns)
    
    def _validate_snv_hgvs(self, line: str) -> SNVInput:
        """Validate HGVS format SNV."""
        parts = line.split(':')
        if len(parts) != 2:
            raise ValueError("HGVS format requires transcript:notation")
        
        transcript_id, hgvs_notation = parts
        
        if transcript_id not in self.transcript_lookup:
            raise ValueError(f"Transcript ID {transcript_id} not found in database")
        
        if hgvs_notation.startswith('c.'):
            if not re.match(r'^c\.\d+[ATCG]>[ATCG]$', hgvs_notation):
                raise ValueError("Invalid HGVS coding format; expected 'c.725G>A'")
            return SNVInput(transcript_id=transcript_id, hgvs_c=hgvs_notation, build="GRCh38")
        elif hgvs_notation.startswith('p.'):
            if not re.match(r'^p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}$', hgvs_notation):
                raise ValueError("Invalid HGVS protein format; expected 'p.Val600Glu'")
            return SNVInput(transcript_id=transcript_id, hgvs_p=hgvs_notation, build="GRCh38")
        else:
            raise ValueError("HGVS notation must start with c. or p.")
    
    def _validate_snv_genomic(self, line: str) -> SNVInput:
        """Validate genomic coordinate format SNV."""
        match = re.match(r'^(chr\w+):(\d+)([ATCG])>([ATCG])$', line, re.IGNORECASE)
        if not match:
            raise ValueError("Invalid genomic format; expected 'chr7:140753336A>T'")
        
        chromosome, position, ref, alt = match.groups()
        return SNVInput(
            chromosome=chromosome,
            position=int(position),
            reference=ref.upper(),
            alternate=alt.upper(),
            build="GRCh38"
        )
    
    def _validate_cnv(self, line: str) -> CNVInput:
        """Validate CNV format."""
        if ':copy_number=' in line:
            match = re.match(r'^(chr\w+):(\d+)-(\d+):copy_number=([\d.]+)$', line, re.IGNORECASE)
            if not match:
                raise ValueError("Invalid CNV format; expected 'chr17:4123-4127:copy_number=3.5'")
            chromosome, start, end, copy_number = match.groups()
            return CNVInput(
                chromosome=chromosome,
                start=int(start),
                end=int(end),
                copy_number=float(copy_number),
                build="GRCh38"
            )
        else:
            raise ValueError("Invalid CNV format. Missing copy_number. Expected: chr:start-end:copy_number=X")
    
    def _validate_fusion(self, line: str) -> FusionInput:
        """Validate fusion format."""
        match = re.match(r'^(NM_\d+\.\d+|ENST\d+):exon(\w+)-(NM_\d+\.\d+|ENST\d+):exon(\w+)$', line)
        if not match:
            raise ValueError("Fusion requires transcript IDs: e.g. 'NM_xxx:exon1-NM_yyy:exon2'")
        
        transcript5p, exon5p, transcript3p, exon3p = match.groups()
        
        if transcript5p not in self.transcript_lookup:
            raise ValueError(f"5' transcript {transcript5p} not found")
        if transcript3p not in self.transcript_lookup:
            raise ValueError(f"3' transcript {transcript3p} not found")
        
        return FusionInput(
            transcript5p=transcript5p,
            junction5p_exon=exon5p,
            transcript3p=transcript3p,
            junction3p_exon=exon3p,
            build="GRCh38"
        )
    
    def _validate_sv(self, line: str) -> SVInput:
        """Validate structural variant format."""
        match = re.match(r'^(DEL|DUP|INV|TRA|COMPLEX):(.+)$', line, re.IGNORECASE)
        if not match:
            raise ValueError("Invalid SV format; expected 'DEL:chr1:100-200'")
        
        svtype, breakend_info = match.groups()
        
        breakends = []
        if ':' in breakend_info:
            parts = breakend_info.split(':')
            if len(parts) >= 2:
                breakends.append({
                    "chr": parts[0],
                    "pos": int(parts[1].split('-')[0]),
                    "strand": "+",
                    "mate_id": "1"
                })
        
        return SVInput(
            breakends=breakends,
            svtype=svtype.upper(),
            build="GRCh38"
        )
    
    def _validate_tumor_marker(self, line: str) -> TumorMarkerInput:
        """Validate tumor marker format."""
        if re.match(r'^MSI-(H|L)$', line, re.IGNORECASE):
            return TumorMarkerInput(msi_status=line.upper())
        elif re.match(r'^MSS$', line, re.IGNORECASE):
            return TumorMarkerInput(msi_status="MSS")
        elif re.match(r'^TMB:([\d.]+)$', line, re.IGNORECASE):
            score = float(line.split(':')[1])
            return TumorMarkerInput(tmb_score=score)
        elif re.match(r'^HRD:([\d.]+)$', line, re.IGNORECASE):
            score = float(line.split(':')[1])
            return TumorMarkerInput(hrd_score=score)
        elif re.match(r'^(\w+)_expression:([\d.]+)$', line, re.IGNORECASE):
            gene, score = line.split('_expression:')
            return TumorMarkerInput(expression_data={gene: float(score)})
        else:
            raise ValueError("Invalid tumor marker format; expected 'MSI-H', 'TMB:10.5', 'HRD:42.0', or 'GENE_expression:2.5'")
    
    def _load_qc_flags(self) -> Dict[str, Dict[str, str]]:
        """Load QC flags from resources file."""
        try:
            qc_file = os.path.join("resources", "qc_flags", "QC_Flags.tsv")
            if os.path.exists(qc_file):
                df = pd.read_csv(qc_file, sep='\t')
                qc_dict = {}
                for _, row in df.iterrows():
                    qc_dict[row['flag_code']] = {
                        'name': row['flag_name'],
                        'description': row['description'],
                        'severity': row['severity'],
                        'action': row['recommended_action']
                    }
                return qc_dict
            else:
                return {
                    'QNS': {'name': 'Quantity Not Sufficient', 'severity': 'HIGH'},
                    'LOW_DEPTH': {'name': 'Low Read Depth', 'severity': 'MEDIUM'},
                    'FAIL_AMP': {'name': 'Failed to Amplify', 'severity': 'HIGH'},
                    'POOR_QUAL': {'name': 'Poor Quality', 'severity': 'MEDIUM'},
                    'CONTAM': {'name': 'Contamination', 'severity': 'HIGH'}
                }
        except Exception as e:
            print(f"Warning: Could not load QC flags: {e}")
            return {}
    
    def _validate_qc_flags(self, qc_flags: QCFlags) -> List[ValidationError]:
        """Validate QC flags against known flag types."""
        errors = []
        
        flag_fields = ['qns', 'low_read_depth', 'failed_to_amplify', 'poor_quality', 'contamination']
        active_flags = []
        
        for field in flag_fields:
            if getattr(qc_flags, field, None) is True:
                flag_code = field.upper().replace('_', '_')
                if field == 'low_read_depth':
                    flag_code = 'LOW_DEPTH'
                elif field == 'failed_to_amplify':
                    flag_code = 'FAIL_AMP'
                elif field == 'poor_quality':
                    flag_code = 'POOR_QUAL'
                elif field == 'contamination':
                    flag_code = 'CONTAM'
                
                active_flags.append(flag_code)
        
        high_severity_flags = [flag for flag in active_flags if self.qc_flags_config.get(flag, {}).get('severity') == 'HIGH']
        if high_severity_flags:
            errors.append(ValidationError(
                line=0,
                input="qc_flags",
                message=f"High severity QC flags detected: {', '.join(high_severity_flags)}. Consider sample reprocessing."
            ))
        
        return errors
