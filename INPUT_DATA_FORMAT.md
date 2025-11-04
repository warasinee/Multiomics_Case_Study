# Input Data: Preparation and Requirements

This document describes the minimal input files and preparation steps required before running the Multiomics analysis pipeline. It focuses on data preparation and requirements (no detailed tool examples). Replace or adapt file names to match your project.

---

## Purpose

Provide clear, minimal requirements and recommended tools for preparing external input files (data that are not produced by the workflow) so they can be used directly by the analysis pipeline.

---

## Required input files (high level)

- Multi-omics data matrices (feature × sample) for each omics layer you will analyse (transcriptomics, proteomics, metabolomics).
- Differential expression results for enrichment analyses (one combined file across strains/experiments).
- Functional annotation mappings linking features to terms/pathways (GO, KEGG, eggNOG, COG, etc.).
- Optional: network node/edge tables for visualization apps (Shiny).

---

## File format requirements (summary)

- Multi-omics matrices:
  - Format: CSV or TSV, features as rows and samples as columns.
  - Column 1: unique feature identifier (e.g., locus tag, gene ID, protein ID).
  - Remaining columns: sample measurements; column headers must match across omics layers and clearly indicate sample/condition.
  - Numeric values only for measurements (use NA or 0 per your downstream tool requirements).

- Differential expression results (single file):
  - Provide a combined long-format table (R data frame or TSV/CSV) with at least these fields: entity_id, species, strain, experiment, logFC, FDR, is_species_core.
  - This file is used for ORA/GSEA across strains and omics layers.

- Functional annotation mappings:
  - Provide a table with columns: entity_id, entity_type (gene/protein/metabolite), fterm_id, fterm_name, fterm_origin.
  - Remove duplicates and ensure IDs match the feature IDs used in your omics matrices and DE results.

- Network node/edge files (optional):
  - Tab-separated files with clear headers for node attributes (name, counts, p.adjust, etc.) and edge attributes (from, to, similarity, counts).

---

## Microbial functional annotation (recommended tools)

For microbial genomes, use dedicated prokaryotic annotation pipelines rather than organism-specific Bioconductor packages. Recommended tools:

- MicrobeAnnotator — https://github.com/cruizperez/MicrobeAnnotator
- Bakta — https://github.com/oschwengers/bakta
- eggNOG-mapper — https://github.com/eggnogdb/eggnog-mapper

Run one or more of these tools on predicted proteins or assemblies and extract per-feature mappings to GO/KEGG/eggNOG/COG terms. Consolidate mappings into the functional annotation table described above.

---

## Differential expression (recommended viewer)

- Degust (visualization, exploration, and export of differential expression results): https://degust.erc.monash.edu/ — Degust can export DE result tables for downstream analyses.

You may still run standard DE tools (DESeq2, edgeR, limma, MSstats) to compute fold-changes and adjusted p-values, then combine results into the required DE file format for downstream enrichment and visualization.

---

## General data preparation tips (brief)

- Use UTF-8 encoding for text files.
- Use "." as decimal separator.
- Keep feature IDs consistent across all files.
- Keep sample names and column order identical across omics layers.
- Remove or clearly mark missing values (NA) as required by downstream tools.
- Document any preprocessing (normalization, transformation) applied to the data; the pipeline expects inputs already prepared appropriately for the chosen methods.

---

## Minimal sample size guidance

- Recommended: at least 6 biological replicates per condition for basic DE testing; more replicates improve power and robustness.

---

## Where to get help

- See the repository README for pipeline-specific instructions and example input file locations.
- For annotation tool usage and installation, consult each tool's GitHub page (links above).
- For DE visualization, consult Degust documentation and website.

---

**Last Updated**: 2025-11-03