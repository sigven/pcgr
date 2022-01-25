# pcgr dev (January 2022)

- Complete restructure of Python and R components.
  - Installation now relies on two separate [conda](https://docs.conda.io/en/latest/)
    packages, `pcgr` (Python component) and `pcgrr` (R component).
    Direct Docker support remains, with the Dockerfile simplified to rely
    exclusively on the installation of the above Conda packages.
- Remove VCF validation step. Feedback from users suggested that Ensembl's
  `vcf-validator` was often too stringent so its use has been deprecated.
  The `--no_vcf_validate` option remains for backwards compatibility.

# pcgr 0.9.2 (June 2021)

- Data bundle updates (CIViC, ClinVar, CancerMine, UniProt KB, PFAM)
- Software upgrades: VEP (104), R v4.1/BioConductor 3.13
- **NEW**: TOML configuration removed - all options to PCGR are now command-line based
- **NEW**: Feed PCGR with a [CPSR report](https://github.com/sigven/cpsr) to view key germline findings in the tumor report
- [CHANGELOG](http://pcgr.readthedocs.io/en/latest/CHANGELOG.html)
- Planned for next release: Support for analysis of RNA fusions


# pcgr 0.9.1 (November 2020)

- Data bundle updates (CIViC, ClinVar, CancerMine, UniProt KB)
- [CHANGELOG](http://pcgr.readthedocs.io/en/latest/CHANGELOG.html)
