# Data Directory

This directory contains example data files for the multi-omics analysis pipeline.

## 📋 Input File Formats

For detailed specifications of all input file formats, see: **[../INPUT_DATA_FORMAT.md](../INPUT_DATA_FORMAT.md)**

## Example Files Provided

### Multi-Omics Data (for PCA & DIABLO)
- `Strep_transcriptome_data.csv` - Gene expression data (5,186 genes × 60 samples)
- `Strep_proteome_data.csv` - Protein abundance data (1,234 proteins × 60 samples)
- `Strep_metabolome_GC_data.csv` - Metabolite concentration data (136 metabolites × 60 samples)

### Differential Expression Results (for ORA/GSEA)
- `All_entity_fold_change_long_format.rds` - DE results across all omics and strains
- `All_funct_annot_clean_noduplicate.rda` - Functional annotations (GO, KEGG, etc.)

### Enrichment Analysis Outputs
- `ORA_GSEA_list_Strep.rds` - Pre-computed ORA and GSEA results
- `Strep_Trans_GSEA_list.rds` - Transcriptomics GSEA results

### Network Files (for visualization)
- `Strep_all_node_trans_GSEA_KEGG.tsv` - Node table for network visualization
- `Strep_all_edge_trans_GSEA_KEGG.tsv` - Edge table for network visualization
- `Strep_all_node_DIABLO_trans_GSEA_KEGG.tsv` - DIABLO network nodes
- `Strep_all_edge_DIABLO_trans_GSEA_KEGG.tsv` - DIABLO network edges

### Annotation Files
- `Strep_gene_anno.tsv` - Gene functional annotations
- `Strep_GO_anno.tsv` - Gene Ontology annotations
- `Strep_KEGG_anno.tsv` - KEGG pathway annotations

### DIABLO & Network Outputs
- `final.plsda.diablo.rds` - Final DIABLO model object
- `myNetwork_conserved_Strep.gml` - Network file for Cytoscape (cutoff=0.9)
- `myNetwork_conserved_Strep_05.gml` - Network file for Cytoscape (cutoff=0.5)
- `my.network.Strep.220425.cys` - Cytoscape session file

## Usage

These files serve as:
1. **Templates** - Examples of proper file formatting
2. **Test Data** - For running the analysis pipeline
3. **Reference** - For validating your own data preparation

## Need Help?

- **Format questions**: See [INPUT_DATA_FORMAT.md](../INPUT_DATA_FORMAT.md)
- **Pipeline questions**: See [README.md](../README.md)
- **Issues**: Open an issue on GitHub
