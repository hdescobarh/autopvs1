# Changelog
All notable changes to this project will be documented in this file.


## 0.1.0 — Unreleased
- Package the fork with uv, a lockfile, and a src layout; require Python 3.11+.
- Move bundled maxentpy and pyhgvs into the private autopvs1._vendor namespace.
- Require AUTOPVS1_CONFIG; resolve reference paths relative to that INI file,
  expanding environment variables and home-directory paths.
- Include vendor resources in distribution artifacts; keep reference datasets external.
- Preserve public classes and defer configuration and reference loading until use.
- Support hg19-only and hg38-only installations with cached assembly bundles,
  shared reference tables, and lazy compatibility aliases in read_data.
- Reject unsupported assembly names consistently and report contextual reference
  errors before launching VEP; failed loads can be retried.
- Add synthetic import, configuration, resource-lifetime, analysis integration,
  and vendor-resource checks.

## 2021-06-30
- Major update: Compatible with hg19/GRCh37 and hg38/GRCh38
- VEP version upgrade to release/104
- Transcript structure update to UCSC ncbiRefSeq 
    * hg19 ncbiRefSeq.txt.gz 20210518
    * hg38 ncbiRefSeq.txt.gz 20210201
- Population allele frequency database
    * hg19: gnomAD r2.1 exomes and genomes
    * hg38: gnomAD r3.0 genomes
- Clinvar version update to archive 2021-06
    * Biologically-relevant transcripts update
    * Functional region and hotspot update
    * Pathogenic sites update
    * Exon with frequent LoF update
- PVS1 levels
    * ClinGen clinical validity update to 20210624
    * Null mouse model update to 20210510


## 2020-07-15
- VEP version upgrade to release/100.
- Transcript structure update to UCSC 20200301.
- Clinvar version update to 20200629.
	* Functional region and hotspot update.
	* Pathogenic sites update.
- PVS1 levels update.
	* ClinGen clinical validity update to 20200630.
	* Null mouse model from IMPC update.


## 2020-03-16
- Support custom transcript.


## 2020-02-08
- Clinvar version update to 20200106.
