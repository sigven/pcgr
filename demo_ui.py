import streamlit as st
import requests
import json
from typing import List, Dict, Any
import re

st.set_page_config(
    page_title="DGG Variant Interpretation Engine",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 DGG Variant Interpretation Engine")
st.markdown("**Clinical variant classification and interpretation demo**")
st.markdown("---")

API_BASE_URL = "http://localhost:8000"
CLASSIFY_ENDPOINT = f"{API_BASE_URL}/v1/classify"

def validate_variant_input(variant_text: str, tumor_type: str) -> Dict[str, Any]:
    """Call the enhanced validation API."""
    variant_lines = [line.strip() for line in variant_text.split('\n') if line.strip()]
    
    try:
        response = requests.post(
            f"{API_BASE_URL}/v1/validate_variants",
            json={
                "variant_lines": variant_lines,
                "oncotree_code": tumor_type
            },
            headers={"Content-Type": "application/json"},
            timeout=30
        )
        
        if response.status_code == 400:
            return response.json()
        
        response.raise_for_status()
        return response.json()
        
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 400:
            return e.response.json()
        raise
    except requests.exceptions.ConnectionError:
        st.error("❌ Could not connect to the DGG Engine API. Make sure the FastAPI server is running on localhost:8000")
        return None
    except Exception as e:
        st.error(f"❌ Validation error: {str(e)}")
        return None


def parse_variant_input(variant_text: str) -> List[Dict[str, Any]]:
    """Parse simple variant text input into structured Variant objects.
    
    Supports formats like:
    - BRAF V600E (creates SNVInput)
    - TP53 amplification (creates CNVInput)
    - BCR-ABL1 fusion (creates SVInput)
    - chr7:140753336A>T (creates SNVInput)
    """
    variants = []
    lines = [line.strip() for line in variant_text.split('\n') if line.strip()]
    
    for line in lines:
        aa_pattern = r'^([A-Z0-9]+)\s+([A-Z]\d+[A-Z])$'
        aa_match = re.match(aa_pattern, line.upper())
        
        if aa_match:
            gene_symbol, aa_change = aa_match.groups()
            variant = {
                "chromosome": "chr1",
                "position": 100000,
                "reference_allele": aa_change[0],
                "alternate_allele": aa_change[-1]
            }
            variants.append(variant)
            continue
        
        cnv_pattern = r'^([A-Z0-9]+)\s+(amplification|loss|deletion|high-level amplification)$'
        cnv_match = re.match(cnv_pattern, line, re.IGNORECASE)
        
        if cnv_match:
            gene_symbol, cnv_type = cnv_match.groups()
            variant = {
                "gene": gene_symbol.upper(),
                "copy_number_change": cnv_type
            }
            variants.append(variant)
            continue
        
        fusion_pattern = r'^([A-Z0-9]+)-([A-Z0-9]+)\s+(fusion|translocation)$'
        fusion_match = re.match(fusion_pattern, line, re.IGNORECASE)
        
        if fusion_match:
            gene1, gene2, sv_type = fusion_match.groups()
            variant = {
                "sv_type": sv_type.lower(),
                "genes": [gene1, gene2]
            }
            variants.append(variant)
            continue
        
        genomic_pattern = r'^chr(\w+):(\d+)([ATCG])>([ATCG])$'
        genomic_match = re.match(genomic_pattern, line, re.IGNORECASE)
        
        if genomic_match:
            chromosome, position, ref, alt = genomic_match.groups()
            variant = {
                "chromosome": f"chr{chromosome}",
                "position": int(position),
                "reference_allele": ref,
                "alternate_allele": alt
            }
            variants.append(variant)
            continue
        
        st.warning(f"Could not parse variant format: {line}. Creating basic SNV.")
        variant = {
            "chromosome": "chr1",
            "position": 100000,
            "reference_allele": "A",
            "alternate_allele": "T"
        }
        variants.append(variant)
    
    return variants

def call_classify_api(variants: List[Dict[str, Any]], tumor_type: str) -> Dict[str, Any]:
    """Make API call to the classify endpoint."""
    request_data = {
        "variants": variants,
        "tumor_type": tumor_type
    }
    
    try:
        response = requests.post(
            CLASSIFY_ENDPOINT,
            json=request_data,
            headers={"Content-Type": "application/json"},
            timeout=30
        )
        response.raise_for_status()
        return response.json()
    except requests.exceptions.ConnectionError:
        st.error("❌ Could not connect to the DGG Engine API. Make sure the FastAPI server is running on localhost:8000")
        st.code("uvicorn api.main:app --reload --host 0.0.0.0 --port 8000")
        return None
    except requests.exceptions.Timeout:
        st.error("⏱️ API request timed out. The classification may be taking longer than expected.")
        return None
    except requests.exceptions.HTTPError as e:
        st.error(f"❌ API request failed with status {e.response.status_code}")
        if e.response.text:
            st.code(e.response.text)
        return None
    except Exception as e:
        st.error(f"❌ Unexpected error: {str(e)}")
        return None

def display_classification_results(results: Dict[str, Any]):
    """Display the classification results in a user-friendly format."""
    if not results:
        return
    
    st.success("✅ Classification completed successfully!")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Variants", results.get("total_variants", 0))
    with col2:
        st.metric("Tumor Type", results.get("tumor_type", "Unknown"))
    with col3:
        pertinent_negatives = results.get("pertinent_negatives", [])
        st.metric("Pertinent Negatives", len(pertinent_negatives))
    
    st.markdown("---")
    
    classified_variants = results.get("classified_variants", [])
    
    if not classified_variants:
        st.warning("No classified variants returned.")
        return
    
    st.subheader("🎯 Classification Results")
    
    for i, variant in enumerate(classified_variants, 1):
        with st.expander(f"Variant {i}: {variant.get('original_variant', {}).get('gene_symbol', 'Unknown')} - Tier {variant.get('tier', 'Unknown')}", expanded=True):
            
            col1, col2 = st.columns([1, 2])
            
            with col1:
                tier = variant.get("tier", "Unknown")
                tier_description = variant.get("tier_description", "No description available")
                
                tier_colors = {
                    "1": "🔴",
                    "2": "🟠", 
                    "3": "🟡",
                    "4": "🔵",
                    "5": "⚪"
                }
                tier_icon = tier_colors.get(str(tier), "❓")
                
                st.markdown(f"### {tier_icon} Tier {tier}")
                st.markdown(f"**Description:** {tier_description}")
            
            with col2:
                orig_variant = variant.get("original_variant", {})
                st.markdown("**Variant Details:**")
                st.markdown(f"- **Gene:** {orig_variant.get('gene_symbol', 'Unknown')}")
                st.markdown(f"- **Type:** {orig_variant.get('variant_type', 'Unknown')}")
                st.markdown(f"- **HGVS Protein:** {orig_variant.get('hgvs_protein', 'N/A')}")
                st.markdown(f"- **HGVS Genomic:** {orig_variant.get('hgvs_genomic', 'N/A')}")
            
            vep_annotation = variant.get("vep_annotation")
            if vep_annotation:
                st.markdown("**🔬 VEP Annotation:**")
                annotation_cols = st.columns(3)
                
                with annotation_cols[0]:
                    st.markdown(f"- **Consequence:** {vep_annotation.get('consequence', 'Unknown')}")
                    st.markdown(f"- **Impact:** {vep_annotation.get('impact', 'Unknown')}")
                
                with annotation_cols[1]:
                    st.markdown(f"- **gnomAD AF:** {vep_annotation.get('gnomad_af', 'N/A')}")
                    st.markdown(f"- **COSMIC ID:** {vep_annotation.get('cosmic_id', 'N/A')}")
                
                with annotation_cols[2]:
                    st.markdown(f"- **dbSNP ID:** {vep_annotation.get('dbsnp_rsid', 'N/A')}")
                    st.markdown(f"- **SIFT:** {vep_annotation.get('sift_prediction', 'N/A')}")
            
            evidence_summary = variant.get("evidence_summary", [])
            if evidence_summary:
                st.markdown("**📋 Evidence Summary:**")
                if isinstance(evidence_summary, list):
                    evidence_text = ", ".join(evidence_summary)
                else:
                    evidence_text = str(evidence_summary)
                st.markdown(f"- {evidence_text}")
    
    if pertinent_negatives:
        st.markdown("---")
        st.subheader("⚠️ Pertinent Negatives")
        for negative in pertinent_negatives:
            st.markdown(f"- {negative}")

st.sidebar.header("🔧 Configuration")

tumor_types = [
    "MEL", "LUAD", "BRCA", "CRC", "PRAD", "BLCA", "HNSC", "KIRC", 
    "LIHC", "THCA", "UCEC", "KIRP", "SARC", "LAML", "PAAD", "GBM",
    "OV", "SKCM", "CESC", "TGCT", "UCS", "CHOL", "THYM", "ACC",
    "MESO", "UVM", "DLBC", "KICH", "READ", "LGG", "PCPG", "STAD"
]

tumor_type = st.sidebar.selectbox(
    "Select Tumor Type",
    options=tumor_types,
    index=0,
    help="OncoTree tumor type code for classification context"
)

st.subheader("📝 Variant Input")
st.markdown("Enter variants in one of these formats:")
st.markdown("- **HGVS Notation:** `NM_004333.4:c.1799T>A` or `NM_004333.4:p.Val600Glu`")
st.markdown("- **Genomic Coordinate:** `chr7:140753336A>T`")
st.markdown("- **CNV:** `chr17:4123-4127:copy_number=3.5`")
st.markdown("- **Fusion:** `NM_001:exon1-NM_002:exon2`")
st.markdown("- **Tumor Markers:** `MSI-H`, `TMB:15.2`, `HRD:42.0`, `BRAF_expression:2.5`")

variant_input = st.text_area(
    "Variants (one per line)",
    value="NM_004333.4:c.1799T>A\nchr7:140753336A>T\nchr17:4123-4127:copy_number=3.5\nMSI-H",
    height=150,
    help="Enter one variant per line using the supported formats above"
)

col1, col2, col3 = st.columns([1, 1, 1])
with col2:
    submit_button = st.button("🚀 Classify Variants", type="primary", use_container_width=True)

if submit_button:
    if not variant_input.strip():
        st.error("❌ Please enter at least one variant to classify.")
    else:
        with st.spinner("🔄 Validating variants..."):
            validation_result = validate_variant_input(variant_input, tumor_type)
            
            if validation_result is None:
                st.stop()
            
            if "errors" in validation_result and validation_result["errors"]:
                st.error("❌ Validation errors found:")
                
                for error in validation_result["errors"]:
                    st.error(f"Line {error['line']}: {error['message']}")
                    st.code(f"Input: {error['input']}")
                
                st.markdown("### 📋 Expected Formats:")
                st.markdown("- **HGVS:** `NM_004333.4:c.1799T>A`")
                st.markdown("- **Genomic:** `chr7:140753336A>T`")
                st.markdown("- **CNV:** `chr17:4123-4127:copy_number=3.5`")
                st.markdown("- **Fusion:** `NM_001:exon1-NM_002:exon2`")
                st.markdown("- **Tumor Markers:** `MSI-H`, `TMB:15.2`")
                
            else:
                st.success("✅ All variants validated successfully!")
                
                summary = validation_result.get("summary", {})
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total Variants", summary.get("total_variants", 0))
                with col2:
                    st.metric("Parsed Successfully", summary.get("parsed", 0))
                with col3:
                    st.metric("Validation Errors", summary.get("errors", 0))
                
                parsed_variants = parse_variant_input(variant_input)
                results = call_classify_api(parsed_variants, tumor_type)
                
                if results:
                    display_classification_results(results)

st.markdown("---")
st.markdown("**DGG Variant Interpretation Engine** - Clinical variant classification demo")
st.markdown("Make sure the FastAPI server is running: `uvicorn api.main:app --reload --host 0.0.0.0 --port 8000`")
