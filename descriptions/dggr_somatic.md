Of course. The provided specification is incredibly detailed but formatted for human comprehension. To optimize it for an agentic AI coding system, we need to translate the prose into a structured, hierarchical, and actionable blueprint.

Here is a reworked specification designed to be parsed and executed by an AI agent. It breaks the project into distinct, implementable components with clear data schemas, function signatures, and logic flows.

***

## Project Specification: DGG Variant Tiering & Canned-Text Engine

**Objective:** Develop a Python-based system to automatically classify genomic variants from solid tumor NGS panels and generate standardized interpretive text for clinical reporting.

**Core Technologies:**
* **Language:** Python 3.10+
* **API Framework:** FastAPI
* **Data Validation:** Pydantic 2.0
* **Core Libraries:** Pandas
* **UI Demo:** Streamlit
* **Containerization:** Docker
* **External Dependency:** ANNOVAR (annotation tool)

---

### **Tier 1: System Architecture & Project Structure**

The system consists of a REST API, a core tiering engine, a knowledge base, and an external annotation tool wrapper.

**Suggested Project Directory Structure:**

```
dgg_rules_engine/
├── api/
│   ├── __init__.py
│   ├── main.py             # FastAPI application
│   └── routes.py           # API routes/endpoints
├── core/
│   ├── __init__.py
│   ├── engine.py           # Main classification logic
│   ├── models.py           # Pydantic data models
│   ├── kb_manager.py       # Knowledge base loader
│   └── pertinent_negatives.py # Logic for negative findings
├── wrappers/
│   ├── __init__.py
│   └── annovar_wrapper.py  # Wrapper for ANNOVAR tool
├── data/
│   ├── knowledge_base/     # Houses the 9 flat-file KBs
│   │   ├── general_gene_text.csv
│   │   └── ... (8 other files)
│   └── builds/             # Output from the ETL process
├── etl/
│   └── build_reportable.py # Script to build variant catalogue
├── tests/
│   ├── __init__.py
│   ├── test_engine.py
│   ├── test_api.py
│   └── ...
├── demo_ui.py                # Streamlit UI application
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

---

### **Tier 2: Component Implementation Details**

#### **Component 1: Data Models (`core/models.py`)**

**Action:** Define Pydantic models to standardize all data structures.

1.  **FeatureType (Enum):**
    * `SNV`, `INDEL`, `CNA`, `FUSION`

2.  **BaseFeature (Pydantic Model):**
    * `feature_type`: `FeatureType`
    * `gene_symbol`: `str`
    * `raw_input`: `dict` (to store the original input)

3.  **Specific Feature Models (inherit from BaseFeature):**
    * **SNVFeature:** `chromosome`: `str`, `position`: `int`, `ref_allele`: `str`, `alt_allele`: `str`
    * **IndelFeature:** `chromosome`: `str`, `start`: `int`, `end`: `int`, `ref_allele`: `str`, `alt_allele`: `str`
    * **CNAFeature:** `copy_number_value`: `Optional[float]`, `qualitative_flag`: `str` (e.g., "Amplified")
    * **FusionFeature:** `gene_5_prime`: `str`, `gene_3_prime`: `str`, `breakpoint_details`: `Optional[str]`

4.  **AnnotationRecord (Pydantic Model):**
    * `gene`: `str`
    * `variant_function`: `str` (e.g., "nonsynonymous SNV")
    * `amino_acid_change`: `Optional[str]`
    * `gnomad_af`: `Optional[float]` (population frequency)
    * `cosmic_id`: `Optional[str]`
    * `clinvar_significance`: `Optional[str]`

5.  **ClassifiedVariant (Pydantic Model):**
    * `feature`: `dict` (the normalized input feature)
    * `annotation`: `AnnotationRecord`
    * `tier`: `str` (e.g., "I", "II", "III", "IV")
    * `interpretation_text`: `str`

6.  **ClassificationResult (Pydantic Model):**
    * `classified_variants`: `List[ClassifiedVariant]`
    * `pertinent_negatives`: `List[str]`
    * `technical_comments`: `List[str]`

#### **Component 2: Knowledge Base (`core/kb_manager.py`)**

**Action:** Implement a `KnowledgeBase` class to load and provide access to the 9 curated data files from `data/knowledge_base/`. All files should be loaded into Pandas DataFrames.

* **`KnowledgeBase` Class:**
    * `__init__(self, kb_path)`: Loads all 9 CSV files into instance attributes (e.g., `self.general_gene_text_df`).
    * Implement getter methods for each data source (e.g., `get_gene_text(gene_symbol)`, `get_variant_summary(variant_id)`).

* **Knowledge Base Schemas (as CSV columns):**
    1.  `general_gene_text.csv`: `gene_symbol`, `gene_name`, `summary_text`, `references`
    2.  `gene_dx_interpretations.csv`: `gene_symbol`, `diagnosis`, `dx_text`, `tier`, `references`
    3.  `variant_summaries.csv`: `variant_id`, `gene_symbol`, `variant_protein`, `oncogenic_category`, `summary_text`, `sources`
    4.  `variant_dx_interpretations.csv`: `variant_id`, `diagnosis`, `tier`, `interpretation_text`, `therapy_association`, `references`
    5.  `incidental_findings.csv`: `gene_symbol`, `variant_id`, `incidental_text`, `report_action`, `references`
    6.  `chr_alterations.csv`: `alteration_id`, `category`, `genes_involved`, `interpretation_text`, `tier`, `references`
    7.  `technical_comments.csv`: `flag_name`, `severity`, `comment_text`, `recommendation`
    8.  `pertinent_negatives.csv`: `diagnosis`, `gene_list`, `statement_template`, `priority`, `references`
    9.  `tmb_msi_comments.csv`: `biomarker`, `category`, `threshold`, `comment_text`, `references`

#### **Component 3: ANNOVAR Wrapper (`wrappers/annovar_wrapper.py`)**

**Action:** Create a function to interface with the ANNOVAR command-line tool. This function should be called from within a Docker container that has ANNOVAR installed.

* **Function Signature:** `run_annovar(features: List[Union[SNVFeature, IndelFeature]]) -> List[AnnotationRecord]`
* **Logic:**
    1.  Accept a list of `SNVFeature` and `IndelFeature` objects.
    2.  Convert the list into ANNOVAR's required input format (AVinput).
    3.  Execute the `table_annovar.pl` script as a subprocess.
    4.  Parse the resulting `multianno.txt` output file.
    5.  For each input variant, create and return a corresponding `AnnotationRecord` Pydantic model.
    6.  Handle errors gracefully if ANNOVAR fails.

#### **Component 4: Core Tiering Engine (`core/engine.py`)**

**Action:** Implement the main classification workflow.

* **Main Function:** `classify_features(features: List[dict], tumor_type: str, kb: KnowledgeBase) -> ClassificationResult`
* **Algorithm:**
    1.  **Normalize & Validate:**
        * Iterate through the raw `features` list.
        * For each feature, determine its type (`SNV`, `INDEL`, etc.) and convert it into the corresponding Pydantic model from `core.models`.
        * Validate required fields. Invalid features should be logged and skipped.
    2.  **Annotate:**
        * Separate SNVs/Indels and pass them to `annovar_wrapper.run_annovar()`.
        * Merge the returned `AnnotationRecord` list back with the feature objects.
    3.  **Tier & Interpret:**
        * Initialize empty lists for `classified_variants` and `technical_comments`.
        * For each annotated feature, call helper functions:
            * `tier = _assign_tier(feature, annotation, kb, tumor_type)`
            * `text = _generate_canned_text(feature, tier, annotation, kb)`
            * Create `ClassifiedVariant` object and append to the list.
    4.  **Generate Pertinent Negatives:**
        * `negatives = pertinent_negatives.generate(detected_genes, tumor_type, kb)`
    5.  **Assemble & Return:**
        * Create and return a `ClassificationResult` object containing the results.

* **Helper Function:** `_assign_tier(...)`
    * **Logic (Decision Tree):**
        1.  **Check Knowledge Base First:** Does the exact variant exist in `variant_dx_interpretations` or `variant_summaries` with a pre-assigned tier? If yes, use it.
        2.  **Apply Rule-Based Logic (AMP/ASCO/CAP style):**
            * **Tier I/II:** Is the variant a known activating mutation in an oncogene? A truncating mutation in a tumor suppressor? A known fusion? Use `gene_level_rules` and annotation data (`variant_function`).
            * **Tier IV:** Is the `gnomad_af > 0.01` (or another threshold)? If yes, assign Tier IV (Benign/Likely Benign).
            * **Tier III:** If none of the above rules apply, default to Tier III (Variant of Uncertain Significance).

* **Helper Function:** `_generate_canned_text(...)`
    * **Logic:**
        1.  **Check for Curated Text:** Look up the variant/gene/tier combination in the `variant_dx_interpretations`, `variant_summaries`, or other relevant KB tables. If a specific text exists, use it.
        2.  **Use Generic Template:** If no specific text is found, generate a standardized text based on the `tier` (e.g., "A variant of uncertain significance was detected in GENE...").

#### **Component 5: REST API (`api/main.py`, `api/routes.py`)**

**Action:** Implement the FastAPI application to expose the engine's functionality.

* **Main Endpoint (Batch Processing):**
    * **Route:** `POST /v1/classify`
    * **Request Body:** JSON object containing `features: List[dict]` and `tumor_type: str`.
    * **Response Body:** A `ClassificationResult` JSON object.
    * **Action:** This endpoint will be the primary entry point, calling `core.engine.classify_features`.

* **Secondary Endpoints (for individual lookups):**
    * `GET /v1/interpretations/gene/{gene_symbol}`
    * `GET /v1/interpretations/variant/{variant_id}`
    * `GET /v1/interpretations/tmb?value={tmb_value}`
    * `GET /v1/interpretations/msi?status={msi_status}`
    * `GET /v1/interpretations/pertinent_negative/{diagnosis}`
    * **Action:** Each of these routes should query the `KnowledgeBase` object to retrieve the corresponding canned text and return it in a structured JSON format.

#### **Component 6: ETL Pipeline (`etl/build_reportable.py`)**

**Action:** Integrate the provided `build_reportable.py` script. This script is responsible for creating the `variant_summaries` and `variant_dx_interpretations` knowledge base files by aggregating data from external sources like OncoKB, CIViC, etc.

* **Purpose:** To create a "Reportable Variant Catalogue".
* **Execution:** This script will be run periodically (e.g., weekly/monthly) via a cron job or manually to keep the internal knowledge base up-to-date with public databases.
* **Output:** Generates TSV/Parquet files that will be used as part of the 9 knowledge base flat files.

#### **Component 7: Containerization (`Dockerfile`, `docker-compose.yml`)**

**Action:** Create Docker configuration to containerize the application and its dependencies.

* **`Dockerfile` for the main application:**
    * Use a `python:3.10-slim` base image.
    * Copy `requirements.txt` and install dependencies.
    * Copy the application code (`api/`, `core/`, etc.).
    * Expose the port for the FastAPI server (e.g., 8000).
    * Use `uvicorn` as the entry point to run the API.

* **`Dockerfile` for ANNOVAR:**
    * Use the Dockerfile provided in "Appendix 1". This container will have ANNOVAR and its databases installed.
    * The `annovar_wrapper` will execute commands inside this container.

* **`docker-compose.yml` (Optional but Recommended):**
    * Define two services: `api_service` and `annovar_service`.
    * This simplifies running the entire stack and manages networking between the services.

---

### **Tier 3: Testing and Demonstration**

#### **Component 8: Testing (`tests/`)**

**Action:** Implement a comprehensive test suite using `pytest`.

* **`test_engine.py`:**
    * Create unit tests for the `_assign_tier` function with mock data to cover every rule.
    * Test the end-to-end `classify_features` function with a known set of inputs and expected outputs.
* **`test_api.py`:**
    * Use FastAPI's `TestClient` to send requests to each endpoint.
    * Assert that the status codes and response bodies are correct.
* **Test Data:** Use small, static CSV files as fixtures to test the `KnowledgeBase` loader, ensuring deterministic test runs.

#### **Component 9: Demo UI (`demo_ui.py`)**

**Action:** Implement the Streamlit demo application.

* **Functionality:**
    1.  Provide input widgets (text boxes, dropdowns) for a user to enter a variant, gene, and diagnosis.
    2.  On submission, use the `requests` library to call the running FastAPI backend (e.g., `POST /v1/classify`).
    3.  Display the returned tier, interpretation text, and pertinent negatives in a clean, readable format.
* **Purpose:** Primarily for user acceptance testing (UAT), demonstration, and knowledge base content review by domain experts.
