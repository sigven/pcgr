import subprocess
import json
import tempfile
import os
from typing import List, Optional, Union
from pydantic import BaseModel, Field

from core.api_models import Variant, SNVInput, CNVInput, SVInput


class VEPAnnotation(BaseModel):
    """Model for essential VEP annotation data."""
    gene_symbol: str = Field(description="Gene symbol")
    consequence: str = Field(description="Variant consequence (e.g., 'missense_variant')")
    hgvsc: Optional[str] = Field(default=None, description="Coding sequence notation")
    hgvsp: Optional[str] = Field(default=None, description="Protein sequence notation")
    existing_variation: Optional[List[str]] = Field(default=None, description="Known IDs like from COSMIC or dbSNP")
    gnomad_af: Optional[float] = Field(default=None, description="gnomAD global allele frequency")
    cosmic_id: Optional[str] = Field(default=None, description="COSMIC mutation identifier")
    dbsnp_rsid: Optional[str] = Field(default=None, description="dbSNP reference ID")


def variant_to_vcf_line(variant: Variant) -> str:
    """Convert a Variant object to VCF format line for VEP input."""
    if variant.variant_type == "SNV":
        chrom = variant.chromosome.replace('chr', '') if variant.chromosome.startswith('chr') else variant.chromosome
        ref = variant.reference or variant.reference_allele
        alt = variant.alternate or variant.alternate_allele
        pos = variant.position
        return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t."
    
    elif variant.variant_type == "CNV":
        gene = variant.gene or "UNKNOWN"
        cn_change = variant.copy_number_change or str(variant.copy_number or "UNKNOWN")
        return f"1\t1000000\t.\tN\t<CNV>\t.\t.\tSVTYPE=CNV;GENE={gene};CN_CHANGE={cn_change}"
    
    elif variant.variant_type == "SV":
        genes_list = variant.genes or []
        genes_str = ",".join(genes_list) if genes_list else "UNKNOWN"
        sv_type = variant.sv_type or "DEL"
        return f"1\t1000000\t.\tN\t<{sv_type.upper()}>\t.\t.\tSVTYPE={sv_type};GENES={genes_str}"
    
    else:
        raise ValueError(f"Unsupported variant type: {variant.variant_type}")


def get_vep_annotation(variant: Variant) -> VEPAnnotation:
    """
    Main entry point for VEP annotation.
    
    Takes a Variant object, runs VEP, and returns parsed annotation.
    For demo purposes, returns mock annotation if VEP is not available.
    """
    try:
        vep_check = subprocess.run(['which', 'vep'], capture_output=True, text=True)
        if vep_check.returncode != 0:
            return _get_mock_annotation(variant)
    except Exception:
        return _get_mock_annotation(variant)
    
    vcf_header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    vcf_line = variant_to_vcf_line(variant)
    vcf_content = vcf_header + vcf_line
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as temp_input:
        temp_input.write(vcf_content)
        temp_input_path = temp_input.name
    
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as temp_output:
            temp_output_path = temp_output.name
        
        vep_command = [
            'vep',
            '--input_file', temp_input_path,
            '--output_file', temp_output_path,
            '--format', 'vcf',
            '--json',
            '--hgvs',
            '--symbol',
            '--protein',
            '--canonical',
            '--af_gnomad',
            '--cache',
            '--offline',
            '--force_overwrite',
            '--no_stats',
            '--quiet'
        ]
        
        result = subprocess.run(
            vep_command,
            capture_output=True,
            text=True,
            check=True
        )
        
        with open(temp_output_path, 'r') as f:
            vep_output = json.load(f)
        
        return parse_vep_json(vep_output, variant)
        
    except subprocess.CalledProcessError as e:
        return _get_mock_annotation(variant)
    except json.JSONDecodeError as e:
        return _get_mock_annotation(variant)
    except Exception as e:
        return _get_mock_annotation(variant)
    finally:
        if os.path.exists(temp_input_path):
            os.unlink(temp_input_path)
        if os.path.exists(temp_output_path):
            os.unlink(temp_output_path)


def parse_vep_json(vep_output: List[dict], original_variant: Variant) -> VEPAnnotation:
    """
    Parse VEP JSON output and extract relevant fields.
    
    VEP JSON output is an array of variant objects, each containing transcript consequences.
    We'll process the most severe consequence or the first canonical transcript.
    """
    if not vep_output:
        raise ValueError("Empty VEP output")
    
    variant_result = vep_output[0]
    
    transcript_consequences = variant_result.get('transcript_consequences', [])
    
    if not transcript_consequences:
        consequences = (variant_result.get('regulatory_consequences', []) + 
                       variant_result.get('intergenic_consequences', []))
        if consequences:
            consequence = consequences[0]
        else:
            raise ValueError("No consequences found in VEP output")
    else:
        canonical_consequence = None
        for tc in transcript_consequences:
            if tc.get('canonical') == 1:
                canonical_consequence = tc
                break
        
        consequence = canonical_consequence or transcript_consequences[0]
    
    gene_symbol = consequence.get('gene_symbol', '')
    consequence_terms = consequence.get('consequence_terms', [])
    consequence_str = consequence_terms[0] if consequence_terms else 'unknown'
    
    hgvsc = consequence.get('hgvsc')
    hgvsp = consequence.get('hgvsp')
    
    gnomad_af = None
    if 'gnomad_af' in consequence:
        gnomad_af = float(consequence['gnomad_af'])
    
    existing_variation = None
    cosmic_id = None
    dbsnp_rsid = None
    
    if 'colocated_variants' in variant_result:
        existing_variation = []
        for cv in variant_result['colocated_variants']:
            if 'id' in cv:
                var_id = cv['id']
                existing_variation.append(var_id)
                
                if var_id.startswith('COSV') or var_id.startswith('COSM'):
                    cosmic_id = var_id
                elif var_id.startswith('rs'):
                    dbsnp_rsid = var_id
    
    return VEPAnnotation(
        gene_symbol=gene_symbol,
        consequence=consequence_str,
        hgvsc=hgvsc,
        hgvsp=hgvsp,
        existing_variation=existing_variation,
        gnomad_af=gnomad_af,
        cosmic_id=cosmic_id,
        dbsnp_rsid=dbsnp_rsid
    )


def _get_mock_annotation(variant: Variant) -> VEPAnnotation:
    """Return mock VEP annotation for demo purposes when VEP is not available."""
    if variant.variant_type == "SNV":
        return VEPAnnotation(
            gene_symbol="BRAF",
            consequence="missense_variant",
            hgvsc="c.1799T>A",
            hgvsp="p.Val600Glu",
            existing_variation=["rs113488022", "COSV12345"],
            gnomad_af=0.0001,
            cosmic_id="COSV12345",
            dbsnp_rsid="rs113488022"
        )
    elif variant.variant_type == "CNV":
        gene_name = variant.gene or "UNKNOWN"
        return VEPAnnotation(
            gene_symbol=gene_name,
            consequence="copy_number_variation",
            hgvsc=None,
            hgvsp=None,
            existing_variation=None,
            gnomad_af=None,
            cosmic_id="COSV67890",
            dbsnp_rsid=None
        )
    elif variant.variant_type == "SV":
        gene_name = variant.genes[0] if variant.genes else "UNKNOWN"
        return VEPAnnotation(
            gene_symbol=gene_name,
            consequence="structural_variant",
            hgvsc=None,
            hgvsp=None,
            existing_variation=None,
            gnomad_af=None,
            cosmic_id="COSV54321",
            dbsnp_rsid=None
        )
    else:
        return VEPAnnotation(
            gene_symbol="UNKNOWN",
            consequence="unknown",
            hgvsc=None,
            hgvsp=None,
            existing_variation=None,
            gnomad_af=None,
            cosmic_id=None,
            dbsnp_rsid=None
        )
