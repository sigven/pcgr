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

## Development Setup

```bash
# Install Poetry if not already installed
curl -sSL https://install.python-poetry.org | python3 -

# Install dependencies
poetry install

# Activate virtual environment
poetry shell
```

## Acknowledgements

This project builds upon the foundational work of the PCGR project. We acknowledge and thank the original PCGR development team for their contributions to the field of cancer genomics interpretation.

Original PCGR Repository: https://github.com/sigven/pcgr

## License

MIT License - see LICENSE file for details.
