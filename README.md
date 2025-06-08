# DGG Variant Interpretation Engine

A modern Python-based system for automated classification of genomic variants and generation of clinical reports.

## Overview

The DGG Variant Interpretation Engine provides:
- 5-tier variant classification based on AMP/ASCO/CAP guidelines
- Explainable outputs with evidence codes
- RESTful API for integration
- Comprehensive knowledge base integration
- Clinical report generation

## Institution

**Michigan Medicine**  
Department of Pathology  
Division of Genomics and Genetics (DGG)

## Author

**Vincent Laufer, MD, PhD**  
Email: vincent.a.laufer@gmail.com

## Architecture

- **API Framework**: FastAPI
- **Data Validation**: Pydantic 2.0
- **Database**: PostgreSQL with SQLModel
- **Annotation**: VEP (Variant Effect Predictor)
- **Containerization**: Docker

## API Endpoints

- `POST /v1/classify_variants`: Main batch processing endpoint
- `POST /v1/validate_variant`: Verify input variants are structurally valid
- `GET /v1/interpretations/gene/{gene_symbol}`: Retrieve gene-level interpretations
- `GET /v1/report/{analysis_id}`: Retrieve generated reports

## Installation

### Local Development Setup

Use the automated setup script:

```bash
# Download and run setup script
curl -sSL https://raw.githubusercontent.com/LauferVA/dgg_rules_somatic/main/local_setup.sh | bash

# Or clone and run manually
git clone https://github.com/LauferVA/dgg_rules_somatic.git
cd dgg_rules_somatic
chmod +x local_setup.sh
./local_setup.sh
```

### Docker Deployment

```bash
# Clone repository
git clone https://github.com/LauferVA/dgg_rules_somatic.git
cd dgg_rules_somatic

# Copy environment template
cp .env.example .env
# Edit .env with your database credentials

# Start services
docker-compose up --build
```

## Metadata Standard

This project uses **PEP 621** format for metadata in `pyproject.toml` with Poetry as the dependency manager. All core project metadata is defined in the `[project]` section following modern Python packaging standards.

## Acknowledgements

This project builds upon the foundational work of the PCGR project. We acknowledge and thank the original PCGR development team for their contributions to the field of cancer genomics interpretation.

Original PCGR Repository: https://github.com/sigven/pcgr

## License

MIT License - see LICENSE file for details.
