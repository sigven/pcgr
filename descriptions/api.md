1. Authentication (JWT, opaque refresh)

Method	Path	Purpose
POST	/api/v1/auth/login	Exchange user credentials for access & refresh tokens
POST	/api/v1/auth/refresh	Issue new access token using valid refresh token
POST	/api/v1/auth/logout	Blacklist presented refresh token

// POST /auth/login – request
{
  "username": "vince",
  "password": "••••••••"
}
// 200 OK – response
{
  "access_token": "eyJhbGc…",
  "access_expires": "2025-05-08T21:10:42Z",
  "refresh_token": "dGhpcyBpcyBsb25nZXI…",
  "refresh_expires": "2025-05-09T09:10:42Z"
}


⸻

2. Sample & Raw-Call Ingestion

Method	Path	Purpose
POST	/api/v1/samples	Register specimen & minimal clinical context
POST	/api/v1/samples/{sample_id}/variants	Bulk-load small-variant calls
POST	/api/v1/samples/{sample_id}/fusions	Bulk-load gene-fusion calls
POST	/api/v1/samples/{sample_id}/cnvs	Bulk-load copy-number events
POST	/api/v1/samples/{sample_id}/qc	Push run/QC metrics
GET	/api/v1/samples/{sample_id}	Retrieve aggregate view

// POST /samples
{
  "sample_id": "S12345",
  "patient_id": "P67890",
  "specimen_type": "FFPE_TUMOR",
  "diagnosis": "Melanoma",
  "reference_genome": "GRCh38"
}
// 201 CREATED
{
  "sample_id": "S12345",
  "created_at": "2025-05-08T16:45:07Z"
}


⸻

3. Rule-Driven Interpretation

Method	Path	Purpose
POST	/api/v1/interpret/variants	Interpret variant array or stored sample variants
POST	/api/v1/interpret/fusions	Interpret fusions
POST	/api/v1/interpret/cnvs	Interpret CNVs
GET	/api/v1/interpret/{sample_id}	Composite interpretation object

// POST /interpret/variants
{
  "sample_id": "S12345",
  "variants": [
    {
      "chrom": "7",
      "pos": 1404513136,
      "ref": "T",
      "alt": "A",
      "gene": "BRAF"
    }
  ]
}
// 200 OK
{
  "sample_id": "S12345",
  "interpretations": [
    {
      "gene": "BRAF",
      "variant": "c.1799T>A (p.Val600Glu)",
      "tier": "1A",
      "commentary": "Pathogenic activating mutation observed in ~50 % of cutaneous melanomas…"
    }
  ]
}


⸻

4. Knowledge-Base Lookup

Method	Path	Purpose
GET	/api/v1/knowledge/reportable-variants?diagnosis=Melanoma	List hotspot variants for disease
GET	/api/v1/knowledge/pertinent-negatives?diagnosis=Melanoma	List loci to mention as negative
GET	/api/v1/knowledge/genes/{symbol}	Curated gene metadata


⸻

5. Boiler-Plate / Comment Snippets

Method	Path	Purpose
POST	/api/v1/comments	Generate or retrieve standardized comment
GET	/api/v1/comments/{comment_id}	Fetch stored snippet

// POST /comments
{
  "sample_id": "S12345",
  "context": "low_coverage",
  "depth": 84
}
// 200 OK
{
  "comment_id": "cm_001aa4",
  "text": "The results are indeterminate in exon 11 due to average read depth < 100×…"
}


⸻

6. Report Assembly & Delivery

Method	Path	Purpose
POST	/api/v1/reports	Trigger full report build (PDF + JSON)
GET	/api/v1/reports/{report_id}	Report metadata + signed URLs
GET	/api/v1/reports/{report_id}/pdf	Direct PDF download

// POST /reports
{
  "sample_id": "S12345",
  "template_id": "solid_tumor_v2"
}
// 202 ACCEPTED
{
  "report_id": "RPT-d5f08c0e",
  "status": "queued",
  "estimated_ready_at": "2025-05-08T16:55:00Z"
}

// GET /reports/RPT-d5f08c0e
{
  "report_id": "RPT-d5f08c0e",
  "sample_id": "S12345",
  "status": "ready",
  "generated_at": "2025-05-08T16:52:13Z",
  "pdf_url": "https://…/RPT-d5f08c0e.pdf?Signature=…",
  "json_url": "https://…/RPT-d5f08c0e.json?Signature=…"
}


⸻

7. Admin & Observability

Method	Path	Purpose
GET	/api/v1/health	Liveness + build SHA
GET	/api/v1/metrics	Prometheus metrics
GET	/api/v1/audit?sample_id=S12345	Immutable audit log


⸻

Field Conventions
	•	Snake-case keys
	•	UTC ISO-8601 timestamps (YYYY-MM-DDThh:mm:ssZ)
	•	ID prefixes (RPT-, cm_, etc.)
	•	Optional values omitted rather than null
	•	Array order not significant unless noted
