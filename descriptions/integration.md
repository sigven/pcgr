Here is the complete, unified project plan. It integrates all instructions, including the detailed branding and attribution rules, into a single document for execution.

***

## **Project Plan: DGG Variant Interpretation Engine**

### **1.0 Project Objective & Core Principles**

The primary objective is to develop the **DGG Variant Interpretation Engine**, a robust, Python-based system for the automated classification of genomic variants and generation of clinical reports. This project involves a substantial rework and modernization of the legacy PCGR software into a new, scalable architecture using Python 3.11+, FastAPI, and PostgreSQL.

Execution of this plan must adhere to the following core principles:

* **Principle 1: Establish a New Identity.** The project and all its outputs will be identified as the **DGG Engine**. The associated institution is **Michigan Medicine**, Department of Pathology, Division is DGG, and the author is **Vincent Laufer, MD, PhD**.

* **Principle 2: Comprehensive Rebranding.** All internal code (variables, comments, docstrings) and all generated outputs (reports, logs, file names) must reflect the new identity. This requires systematically scrubbing all references to the legacy tool ("PCGR"), its original authors, and associated organizations.

* **Principle 3: Preserve Scientific Evidence.** The rebranding process **must not** alter formal scientific citations within the knowledge base. The integrity of the underlying scientific evidence is paramount. Use regular expressions to identify and protect bibliographic references (e.g., `...et al.,...`, `(PMID: ...)`), ensuring they remain untouched during the text-cleaning process.

* **Principle 4: Maintain Ethical Attribution.** Proper credit to PCGR project as a foundational work must be provided in the `README.md` file and the LICENSE file as well as any resulting publications. Aside from this, scrub all reference from code base and deliverables the code makes.

---

### **2.0 Phased Implementation Plan**

#### **Phase 1: Foundation & Backend Modernization**

*Objective: Establish the new technical architecture, scrub the initial codebase, and build a rebranded, citation-aware database.*

1.  **Initial Codebase Setup & Rebranding:**
    * Fork the legacy PCGR repository.
    * Perform a project-wide search-and-replace to align with **Principle 2**.
        * Scrub legacy terms: `PCGR`, `DGG Variant Interpretation Engine`, `sigven`.
        * Scrub contributor names from code comments and docstrings.
    * Update all file headers to reflect the new authorship (`Vincent`, `Michigan`).
    * Create the `README.md` file and add the external attribution statement as required by **Principle 4**.

2.  **Technology Stack & Database Schema:**
    * Establish the project structure using the target specification.
    * Define the complete database schema in PostgreSQL using `SQLModel`. The schema must accommodate all data from the nine knowledge base files and biomarker interpretations.
    * Initialize `Alembic` for database version control and create the initial migration.

3.  **Knowledge Base Ingestion & Sanitization:**
    * Develop the ETL script (`etl/build_reportable.py`).
    * This script's core function is to parse the source knowledge base files and populate the PostgreSQL database.
    * Integrate a citation-aware sanitizer function that, for every text field, first identifies and protects formal citations (**Principle 3**) before replacing any remaining legacy tool attributions (**Principle 2**).

---

#### **Phase 2: Core Annotation & Tiering Migration**

*Objective: Re-implement the core variant annotation, classification, and interpretation logic within the new framework.*

1.  **Implement Annotation Wrapper:**
    * Create the `wrappers/annovar_wrapper.py` to interface with the ANNOVAR command-line tool for variant annotation. The wrapper will convert Pydantic models to AVinput and parse ANNOVAR's output back into `AnnotationRecord` models.

2.  **Implement 5-Tier Classification Engine:**
    * In `core/engine.py`, develop the `_assign_tier` helper function to implement the 5-tier classification system based on AMP/ASCO/CAP guidelines.
    * This function will execute queries against the sanitized PostgreSQL knowledge base to determine the appropriate tier based on variant and gene evidence.

3.  **Implement Canned-Text Generation:**
    * Develop the `_generate_canned_text` helper function.
    * This function will perform cascading queries against the database to fetch the most specific interpretation text available (variant-specific, gene-specific, etc.).
    * If no curated text is found, it will generate standardized text based on the variant's tier and gene, using the new **DGG Engine** identity.

---

#### **Phase 3: Integration of Advanced Biomarker Analysis**

*Objective: Replace the legacy R-based biomarker analyses with Python-native solutions.*

1.  **Tumor Mutational Burden (TMB):**
    * Implement a Python function in `core/biomarkers.py` to calculate TMB (non-synonymous mutations/Mb). The coding region size must be a configurable parameter.
    * The TMB score will be used to query the database for the appropriate interpretive text (e.g., TMB-High).

2.  **Microsatellite Instability (MSI):**
    * Implement a wrapper for a standard Python-based MSI prediction tool (e.g., MSIsensor-pro).
    * The wrapper will execute the tool and parse its output to determine MSI status (e.g., MSI-High, MSS).

3.  **Mutational Signatures:**
    * Integrate the `SigProfilerExtractor` library to perform mutational signature analysis against the COSMIC database.
    * Parse the results to identify contributing signatures and retrieve their associated interpretations from the database.

4.  **Update Final Data Model:**
    * Expand the final `ClassificationResult` Pydantic model to include dedicated objects for the TMB, MSI, and Mutational Signature analysis results.

---

#### **Phase 4: API Finalization, Reporting & Deployment**

*Objective: Expose all functionality through a robust API, generate fully rebranded reports, and containerize the application for deployment.*

1.  **Develop Reporting Solution:**
    * Create a professional HTML report template using **Jinja2** (`report_template.html`).
    * The template's design, headers, and footers must be unique and branded exclusively with the **DGG Engine**, **Michigan**, and **Vincent** identity.
    * The template will be populated with data from the `ClassificationResult` object, which includes the citation-preserving text from the database.

2.  **Finalize API and Testing:**
    * Refine the `POST /v1/classify` endpoint and create a new `POST /v1/generate_report` endpoint that returns the final HTML report.
    * Expand the `pytest` suite to include tests for biomarker calculations and, critically, to assert that no legacy branding appears in any API output while formal citations are preserved.

3.  **Finalize Containerization:**
    * Update the `Dockerfile` and `docker-compose.yml` for the complete application stack (API, PostgreSQL, ANNOVAR).
    * Build and tag the final Docker images using a new, non-legacy naming scheme (e.g., `dgg-engine/api`).
    * Document the final, rebranded CLI command for end-users.
