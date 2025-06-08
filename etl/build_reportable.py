import pandas as pd
from sqlmodel import Session
from core.kb_models import (
    GeneGeneralComment, VariantDxInterpretation, PertinentNegative,
    VariantGeneralComment, DrugAssociation, OncoTreeMapping
)
from api.database import engine
from datetime import datetime
import os


def load_gene_general_comments():
    """Load general gene comments from TSV."""
    file_path = 'resources/gene_text/general_gene_comments.tsv'
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found, skipping...")
        return
    
    df = pd.read_csv(file_path, sep='\t')
    with Session(engine) as session:
        for _, row in df.iterrows():
            if pd.isna(row['Summary']) or not str(row['Summary']).strip():
                continue
            
            comment = GeneGeneralComment(
                gene_id=str(row['GeneID']),
                gene_symbol=row['Symbol'],
                comment=str(row['Summary']).strip(),
                created_by='etl_pipeline'
            )
            session.add(comment)
        session.commit()
    print(f"Loaded {len(df)} gene general comments")


def load_variant_dx_interpretations():
    """Load variant disease interpretations from TSV."""
    file_path = 'resources/variant_interps/variant_dx_interpretations.tsv'
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found, skipping...")
        return
    
    df = pd.read_csv(file_path, sep='\t')
    with Session(engine) as session:
        for _, row in df.iterrows():
            if pd.isna(row['molecular_profile']) or pd.isna(row['disease']) or pd.isna(row['evidence_statement']):
                continue
                
            interp = VariantDxInterpretation(
                variant_id=str(row['molecular_profile']).strip(),
                molecular_profile_id=str(row.get('molecular_profile_id', '')) if pd.notna(row.get('molecular_profile_id')) else None,
                disease=str(row['disease']).strip(),
                doid=str(row.get('doid', '')) if pd.notna(row.get('doid')) else None,
                interpretation=str(row['evidence_statement']).strip(),
                evidence_type=str(row.get('evidence_type', '')) if pd.notna(row.get('evidence_type')) else None,
                evidence_direction=str(row.get('evidence_direction', '')) if pd.notna(row.get('evidence_direction')) else None,
                evidence_level=str(row.get('evidence_level', '')) if pd.notna(row.get('evidence_level')) else None,
                significance=str(row.get('significance', '')) if pd.notna(row.get('significance')) else None,
                references=str(row.get('citation', '')) if pd.notna(row.get('citation')) else '',
                citation_id=str(row.get('citation_id', '')) if pd.notna(row.get('citation_id')) else None,
                citation=str(row.get('citation', '')) if pd.notna(row.get('citation')) else None,
                nct_ids=str(row.get('nct_ids', '')) if pd.notna(row.get('nct_ids')) else None,
                evidence_status=str(row.get('evidence_status', '')) if pd.notna(row.get('evidence_status')) else None,
                civic_url=str(row.get('evidence_civic_url', '')) if pd.notna(row.get('evidence_civic_url')) else None,
                created_by='etl_pipeline'
            )
            session.add(interp)
        session.commit()
    print(f"Loaded {len(df)} variant interpretations")


def load_pertinent_negatives():
    """Load pertinent negatives from TSV."""
    file_path = 'resources/pertinent_negatives/Pertinent_Negatives.tsv'
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found, skipping...")
        return
    
    df = pd.read_csv(file_path, sep='\t')
    with Session(engine) as session:
        for _, row in df.iterrows():
            if pd.isna(row['Pertinent Negatives']) or not row['Pertinent Negatives'].strip():
                continue
            
            genes_text = str(row['Pertinent Negatives']).strip()
            genes = [g.strip() for g in genes_text.split(',') if g.strip()]
            
            for gene in genes:
                gene_clean = gene.replace(' mutation', '').replace(' rearrangement', '').replace(' copy number gain', '').replace(' amplification', '').replace(' copy number variant', '').replace(' copy number alteration', '').replace(' exon 14 skipping mutation', '')
                
                neg = PertinentNegative(
                    gene_symbol=gene_clean,
                    disease=row['Tissue'],
                    comment=str(row.get('Notes', '')),
                    pertinent_negatives_text=genes_text,
                    negative_type=row.get('Type'),
                    created_by='etl_pipeline'
                )
                session.add(neg)
        session.commit()
    print(f"Loaded pertinent negatives from {len(df)} rows")


def load_variant_general_comments():
    """Load general variant comments from TSV."""
    file_path = 'resources/variant_interps/general_variant_comments.tsv'
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found, skipping...")
        return
    
    df = pd.read_csv(file_path, sep='\t')
    with Session(engine) as session:
        for _, row in df.iterrows():
            if pd.isna(row['blurb']) or not str(row['blurb']).strip():
                continue
            if pd.isna(row['gene']) or pd.isna(row['variant']):
                continue
                
            comment = VariantGeneralComment(
                variant_id=str(row['variant_id']) if pd.notna(row.get('variant_id')) else '',
                gene=str(row['gene']).strip(),
                variant=str(row['variant']).strip(),
                comment=str(row['blurb']).strip(),
                chromosome=str(row.get('chromosome', '')) if pd.notna(row.get('chromosome')) else None,
                start=int(row['start']) if pd.notna(row.get('start')) else None,
                stop=int(row['stop']) if pd.notna(row.get('stop')) else None,
                variant_types=str(row.get('variant_types', '')) if pd.notna(row.get('variant_types')) else None,
                hgvs_descriptions=str(row.get('hgvs_descriptions', '')) if pd.notna(row.get('hgvs_descriptions')) else None,
                created_by='etl_pipeline'
            )
            session.add(comment)
        session.commit()
    print(f"Loaded {len(df)} variant general comments")


def load_gene_dx_interpretations():
    """Load gene disease interpretations from TSV."""
    file_path = 'resources/gene_text/gene_dx_interpretation.tsv'
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found, skipping...")
        return
    
    df = pd.read_csv(file_path, sep='\t')
    with Session(engine) as session:
        for _, row in df.iterrows():
            if pd.isna(row['gene_symbol']):
                continue
                
            civic_desc = str(row.get('civic_desc', '')) if pd.notna(row.get('civic_desc')) else ''
            oncokb_desc = str(row.get('oncokb_desc', '')) if pd.notna(row.get('oncokb_desc')) else ''
            
            combined_comment = civic_desc
            if civic_desc and oncokb_desc:
                combined_comment += '\n\n' + oncokb_desc
            elif oncokb_desc:
                combined_comment = oncokb_desc
            
            if not combined_comment.strip():
                continue
                
            comment = GeneGeneralComment(
                gene_id=str(row.get('gene_id', '')) if pd.notna(row.get('gene_id')) else '',
                gene_symbol=str(row['gene_symbol']).strip(),
                comment=combined_comment.strip(),
                created_by='etl_pipeline'
            )
            session.add(comment)
        session.commit()
    print(f"Loaded {len(df)} gene disease interpretations")


def load_biomarkers():
    """Load biomarker data from TSV."""
    file_path = 'resources/Biomarkers/Biomarkers.tsv'
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found, skipping...")
        return
    
    df = pd.read_csv(file_path, sep='\t')
    with Session(engine) as session:
        for _, row in df.iterrows():
            drug_assoc = DrugAssociation(
                variant_id=row['Biomarker'],
                disease='ANY',
                drug_name=row['Biomarker'],
                response_type=row['Status'],
                references=str(row.get('Citations', '')),
                evidence_level='B',
                created_by='etl_pipeline'
            )
            session.add(drug_assoc)
        session.commit()
    print(f"Loaded {len(df)} biomarkers")


def main():
    """Run complete ETL pipeline."""
    print("Starting ETL pipeline...")
    
    load_gene_general_comments()
    load_gene_dx_interpretations()
    load_variant_dx_interpretations()
    load_variant_general_comments()
    load_pertinent_negatives()
    load_biomarkers()
    
    print("ETL pipeline complete")


if __name__ == "__main__":
    main()
