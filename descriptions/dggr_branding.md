Of course. That is a critical and nuanced requirement. The goal is to create a new tool identity while respecting and preserving the underlying scientific citations that form the knowledge base. The plan will be updated to include a more sophisticated, targeted scrubbing process.

Here is the updated project plan with the new instructions integrated. The modifications are highlighted in bold.

***

## Project Plan: Unifying and Rebranding PCGR into the DGG Engine (v3)

**Executive Summary:** This project will fork and substantially refactor the legacy PCGR tool into a modern, maintainable Python application, hereafter named the "DGG Engine." The final product will be a FastAPI-based service with a PostgreSQL backend, capable of detailed variant annotation, a 5-tier classification, and advanced biomarker analysis.

A core principle is a complete rebranding of all code and user-facing outputs. This involves a two-pronged approach:
1.  **Code & Tool Attribution Scrubbing:** All mentions of PCGR, its authors (including original project leads and contributors with European or other national origins), and associated institutions will be scrubbed from the application's internal code and self-referential text. These will be replaced with the new identity ("Institution: Michigan", "Author: Vincent").
2.  **Preservation of Scientific Citations:** **Crucially, all formal bibliographic citations to scientific publications within the knowledge base content (e.g., "Smith et al., Nature 2022", "(PMID: 123456)") will be explicitly preserved.** The integrity of the scientific evidence must not be altered.

Proper attribution to the original PCGR work will be maintained externally in the GitHub repository's `README.md`, `LICENSE`, and all academic publications.

---

### **Phase 0: Project Setup & Foundational Rebranding**

*Objective: To systematically scrub the old identity and contributor names from the initial codebase and establish the new project's identity before development begins.*

1.  **Codebase Forking and Renaming:** (No change)
    * Fork, rename repository to `dgg_rules_engine`, rename executables (`dgg_engine.py`), and rename configuration files (`dgg_config.yaml`).

2.  **Automated Code & Asset Scrubbing:** (No change)
    * Project-wide search-and-replace for: `PCGR` -> `DGG Engine`, `Personal Cancer Genome Reporter` -> `DGG Variant Interpretation Engine`, `sigven` -> `michigan_genomics`.

3.  **Targeted Contributor Name Scrubbing:**
    * **Action: Compile a list of contributor names.** Review the original PCGR repository's `AUTHORS` file, `README`, commit history, and documentation to create a list of primary author and contributor names (e.g., Nakken, and other relevant European and international names).
    * **Action: Perform a targeted search-and-replace.** For each name on the list, search the codebase (primarily comments, docstrings, and plain text files).
    * **Rule:** Replace the name **only when it refers to the authorship of the code or the tool itself.** Do not replace names that appear in what might be a citation. When in doubt, leave the name for the more detailed review in Phase 1.

4.  **Update Author and Copyright Information:** (No change)
    * Modify all file headers and docstrings to change the author and copyright information to `Author: Vincent`, `Institution: Michigan`.

5.  **Project Governance & External Attribution:** (No change)
    * Create a `README.md` with a prominent "Acknowledgements" section crediting the original PCGR project and linking to its source.

---

### **Phase 1: Foundation & Backend Modernization**

*Objective: Establish the new Python architecture and database, ensuring all ingested knowledge base content is sanitized according to the nuanced rebranding rules.*

1.  **Technology Stack & Project Scaffolding:** (No change)
2.  **Database Schema & Migration Setup:** (No change)

3.  **ETL for Knowledge Base Ingestion & Rebranding:**
    * Enhance the `etl/build_reportable.py` script with a sophisticated cleaning function.
    * **Action: Implement a citation-aware text sanitizer.** This function will process all text fields (e.g., `summary_text`, `interpretation_text`) before they are inserted into the database.
    * **Logic:**
        * **1. Preserve Formal Citations:** Use regular expressions to identify and protect scientific citations. Any text matching these patterns **will not be altered**.
            * Pattern examples:
                * `\w+\s+et al\.,\s+\w+\s+\d{4}` (e.g., "Nakken et al., Cell 2021")
                * `\(PMID:\s*\d+\)` (e.g., "(PMID: 3456789)")
                * `\(First author \d{4}\)` (e.g., "(Smith 2022)")
        * **2. Replace Tool Attributions:** After protecting citations, run replacements for tool-specific references.
            * Example: `"The PCGR framework suggests..."` becomes `"The DGG Engine framework suggests..."`.
        * **3. Review Informal Mentions:** For any contributor names found that were not part of a formal citation, replace them.
            * Example: `"An interpretation by Nakken..."` (if it exists outside a formal citation) becomes `"An interpretation by the DGG Engine..."`.

4.  **API and KB Manager Implementation:** The `core/kb_manager.py` will now connect to the fully sanitized and citation-preserving PostgreSQL database.

---

### **Phase 2: Core Annotation & Tiering Migration**

*Objective: Re-implement the core variant processing logic on the clean codebase.*

*(This phase proceeds as planned, building upon the sanitized foundation. All implementation will inherently follow the new branding.)*

---

### **Phase 3: Integration of Advanced Biomarker Analysis**

*Objective: Port the R-based biomarker functionalities to Python, ensuring all outputs are free of legacy branding.*

*(This phase proceeds as planned. All generated interpretation text will be pulled from the sanitized database or created using new templates.)*

---

### **Phase 4: API Finalization, Reporting & Deployment**

*Objective: Solidify the API, create a fully rebranded report, and package the application for deployment under its new name.*

1.  **Solidify the API:** (No change)

2.  **Develop a Rebranded Reporting Solution:**
    * **Goal:** Create a report that is functionally equivalent but visually and textually distinct, while respecting scientific data integrity.
    * **Action:** Use **Jinja2** to populate a newly designed `report_template.html`.
    * **Branding:** The report's header/footer and overall design will reflect the "DGG Engine," "Michigan," and "Vincent" identity.
    * **Content Integrity:** The template will be populated with data directly from the `ClassificationResult` object. Since the underlying interpretation text in the database has already been sanitized to preserve citations, the report will correctly display scientific evidence (e.g., "Pathogenic (Smith et al., NEJM 2020)") while attributing the report itself to the DGG Engine.

3.  **Comprehensive Testing:**
    * Expand the `pytest` suite.
    * **Action: Add specific tests to validate the sanitization logic.** Create test cases with strings containing both tool attributions and formal citations, and assert that only the tool attributions are changed.

4.  **Finalize Containerization & Deployment:** (No change)
    * Use new image names (`ghcr.io/your-org/dgg-engine:latest`).
    * Document and use the fully rebranded CLI command for execution.
