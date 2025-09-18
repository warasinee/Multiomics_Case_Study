Preparing node and edge data for network analysis
================

- [Load R packages](#load-r-packages)
- [Import ORA and GSEA results for all strains and all
  omics](#import-ora-and-gsea-results-for-all-strains-and-all-omics)
- [Collect node and edge data](#collect-node-and-edge-data)
- [Save](#save)

## Load R packages

    library(DOSE)

    ## 

    ## DOSE v3.24.2  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    ## 
    ## If you use DOSE in published research, please cite:
    ## Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015, 31(4):608-609

    library(clusterProfiler)

    ## Registered S3 methods overwritten by 'treeio':
    ##   method              from    
    ##   MRCA.phylo          tidytree
    ##   MRCA.treedata       tidytree
    ##   Nnode.treedata      tidytree
    ##   Ntip.treedata       tidytree
    ##   ancestor.phylo      tidytree
    ##   ancestor.treedata   tidytree
    ##   child.phylo         tidytree
    ##   child.treedata      tidytree
    ##   full_join.phylo     tidytree
    ##   full_join.treedata  tidytree
    ##   groupClade.phylo    tidytree
    ##   groupClade.treedata tidytree
    ##   groupOTU.phylo      tidytree
    ##   groupOTU.treedata   tidytree
    ##   is.rooted.treedata  tidytree
    ##   nodeid.phylo        tidytree
    ##   nodeid.treedata     tidytree
    ##   nodelab.phylo       tidytree
    ##   nodelab.treedata    tidytree
    ##   offspring.phylo     tidytree
    ##   offspring.treedata  tidytree
    ##   parent.phylo        tidytree
    ##   parent.treedata     tidytree
    ##   root.treedata       tidytree
    ##   rootnode.phylo      tidytree
    ##   sibling.phylo       tidytree

    ## clusterProfiler v4.6.2  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    library(RCy3)
    library(igraph)

    ## 
    ## Attaching package: 'igraph'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     simplify

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

    library(tidyverse)

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%--%()      masks igraph::%--%()
    ## ✖ dplyr::as_data_frame() masks tibble::as_data_frame(), igraph::as_data_frame()
    ## ✖ purrr::compose()       masks igraph::compose()
    ## ✖ tidyr::crossing()      masks igraph::crossing()
    ## ✖ dplyr::filter()        masks clusterProfiler::filter(), stats::filter()
    ## ✖ dplyr::lag()           masks stats::lag()
    ## ✖ purrr::simplify()      masks igraph::simplify(), clusterProfiler::simplify()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

## Import ORA and GSEA results for all strains and all omics

    ORA_GSEA_list <- readRDS("/Users/wmujchariyak/Desktop/ORA_GSEA_list_Strep.rds") 

    # Import ORA results for all strains and all omics
    Trans.list.ORA <- ORA_GSEA_list[["Transcriptomics RG"]]$ORA 
    Prot.list.ORA <- ORA_GSEA_list[["Proteomics MS1 DDA"]]$ORA 
    Met.GC.list.ORA <- ORA_GSEA_list[["Metabolomics GC-MS"]]$ORA  

    # Import GSEA results for all strains and all omics
    Trans.list.GSEA <- ORA_GSEA_list[["Transcriptomics RG"]]$GSEA 
    Prot.list.GSEA <- ORA_GSEA_list[["Proteomics MS1 DDA"]]$GSEA
    Met.GC.list.GSEA <- ORA_GSEA_list[["Metabolomics GC-MS"]]$GSEA

## Collect node and edge data

    # Prepare Data for network analysis (e.g., node and edge)
    # Example for only Transcriptomics (ORA - KEGG)

    # 1. Prepare data
    # 1.1 Subset ORA - KEGG AND GO (Example for ORA)
    Trans.5448.ORA.GO.sig <- clusterProfiler::filter(Trans.list.ORA[["5448"]], p.adjust < 0.05, grepl("^GO", ID)) 
    Trans.HKU419.ORA.GO.sig <- clusterProfiler::filter(Trans.list.ORA[["HKU419"]], p.adjust < 0.05, grepl("^GO", ID)) 
    Trans.PS003.ORA.GO.sig <- clusterProfiler::filter(Trans.list.ORA[["PS003"]], p.adjust < 0.05, grepl("^GO", ID)) 
    Trans.PS006.ORA.GO.sig <- clusterProfiler::filter(Trans.list.ORA[["PS006"]], p.adjust < 0.05, grepl("^GO", ID)) 
    Trans.SP444.ORA.GO.sig <- clusterProfiler::filter(Trans.list.ORA[["SP444"]], p.adjust < 0.05, grepl("^GO", ID)) 

    Trans.5448.ORA.KEGG.sig <- clusterProfiler::filter(Trans.list.ORA[["5448"]], p.adjust < 0.05, grepl("^map", ID)) 
    Trans.HKU419.ORA.KEGG.sig <- clusterProfiler::filter(Trans.list.ORA[["HKU419"]], p.adjust < 0.05, grepl("^map", ID)) 
    Trans.PS003.ORA.KEGG.sig <- clusterProfiler::filter(Trans.list.ORA[["PS003"]], p.adjust < 0.05, grepl("^map", ID)) 
    Trans.PS006.ORA.KEGG.sig <- clusterProfiler::filter(Trans.list.ORA[["PS006"]], p.adjust < 0.05, grepl("^map", ID)) 
    Trans.SP444.ORA.KEGG.sig <- clusterProfiler::filter(Trans.list.ORA[["SP444"]], p.adjust < 0.05, grepl("^map", ID)) 

    # 1.2 Subset GSEA - KEGG AND GO (Example for GSEA)
    Trans.5448.GSEA.GO.sig <- clusterProfiler::filter(Trans.list.GSEA[["5448"]], p.adjust < 0.05, grepl("^GO", ID))
    Trans.HKU419.GSEA.GO.sig <- clusterProfiler::filter(Trans.list.GSEA[["HKU419"]], p.adjust < 0.05, grepl("^GO", ID))
    Trans.PS003.GSEA.GO.sig <- clusterProfiler::filter(Trans.list.GSEA[["PS003"]], p.adjust < 0.05, grepl("^GO", ID))
    Trans.PS006.GSEA.GO.sig <- clusterProfiler::filter(Trans.list.GSEA[["PS006"]], p.adjust < 0.05, grepl("^GO", ID))
    Trans.SP444.GSEA.GO.sig <- clusterProfiler::filter(Trans.list.GSEA[["SP444"]], p.adjust < 0.05, grepl("^GO", ID))

    Trans.5448.GSEA.KEGG.sig <- clusterProfiler::filter(Trans.list.GSEA[["5448"]], p.adjust < 0.05, grepl("^map", ID)) 
    Trans.HKU419.GSEA.KEGG.sig <- clusterProfiler::filter(Trans.list.GSEA[["HKU419"]], p.adjust < 0.05, grepl("^map", ID)) 
    Trans.PS003.GSEA.KEGG.sig <- clusterProfiler::filter(Trans.list.GSEA[["PS003"]], p.adjust < 0.05, grepl("^map", ID)) 
    Trans.PS006.GSEA.KEGG.sig <- clusterProfiler::filter(Trans.list.GSEA[["PS006"]], p.adjust < 0.05, grepl("^map", ID)) 
    Trans.SP444.GSEA.KEGG.sig <- clusterProfiler::filter(Trans.list.GSEA[["SP444"]], p.adjust < 0.05, grepl("^map", ID)) 

    # 2. Extract node and edge data from Trans.GSEA.KEGG 
          # Check class of obj. - For example, class(Trans.5448.GSEA.KEGG.sig) -> an enrichResult object 

    # 2.1 Calculate pairwise similarity (of the enriched terms using Jaccard’s similarity index)
          # Check added similarity matrix in the termsim slot of enrichment result 
          # For example, View(Trans.5448.GSEA.KEGG.sig_2@termsim)
    Trans.5448.GSEA.KEGG.sig_2 <- enrichplot::pairwise_termsim(Trans.5448.GSEA.KEGG.sig, method = "JC")  
    Trans.HKU419.GSEA.KEGG.sig_2 <- enrichplot::pairwise_termsim(Trans.HKU419.GSEA.KEGG.sig, method = "JC") 
    Trans.PS003.GSEA.KEGG.sig_2 <- enrichplot::pairwise_termsim(Trans.PS003.GSEA.KEGG.sig, method = "JC")  
    Trans.PS006.GSEA.KEGG.sig_2 <- enrichplot::pairwise_termsim(Trans.PS006.GSEA.KEGG.sig, method = "JC") 
    Trans.SP444.GSEA.KEGG.sig_2 <- enrichplot::pairwise_termsim(Trans.SP444.GSEA.KEGG.sig, method = "JC") 

    # 2.2 Create Enrichment map 
    emap_Trans.5448.GSEA.KEGG.sig <- enrichplot::emapplot(Trans.5448.GSEA.KEGG.sig_2)
    emap_Trans.HKU419.GSEA.KEGG.sig <- enrichplot::emapplot(Trans.HKU419.GSEA.KEGG.sig_2)
    emap_Trans.PS003.GSEA.KEGG.sig <- enrichplot::emapplot(Trans.PS003.GSEA.KEGG.sig_2)
    emap_Trans.PS006.GSEA.KEGG.sig <- enrichplot::emapplot(Trans.PS006.GSEA.KEGG.sig_2)
    emap_Trans.SP444.GSEA.KEGG.sig <- enrichplot::emapplot(Trans.SP444.GSEA.KEGG.sig_2)

    # 2.3 Function to convert enrichment map to dataframe
    emap_to_network_df <- function(emap, what = "edges"){
      if (what == "vertices"){
        df <- emap$data %>% 
          attributes %>% # provide data attributes in list
          .$graph %>% # access an attribute named graph (which is in igraph format containing vertices (node) and edges) 
          igraph::as_data_frame(what = "vertices") %>% # convert igraph nodes (vertices) to dataframe
          rename(gene_set_size = "size") %>%
          rename(p.adjust = "color") %>%
          select(name, gene_set_size, p.adjust)
        return(df)
      } else if (what == "edges"){
        df <- emap$data %>%
          attributes %>%
          .$graph %>%
          igraph::as_data_frame(what = "edges") %>%
          rename(similarity = "weight") %>%
          select(from, to, similarity)
        return(df)
      } else {
        return(print("error what should be set to edges or nodes"))
      }
    }

    # 2.4 Concatenate edges df and node df of the 2 strain
    all_edge_trans_GSEA_KEGG <- rbind(emap_to_network_df(emap_Trans.5448.GSEA.KEGG.sig),
                                     emap_to_network_df(emap_Trans.HKU419.GSEA.KEGG.sig),
                                     emap_to_network_df(emap_Trans.PS003.GSEA.KEGG.sig),
                                     emap_to_network_df(emap_Trans.PS006.GSEA.KEGG.sig),
                                     emap_to_network_df(emap_Trans.SP444.GSEA.KEGG.sig)) %>%
      dplyr::rowwise() %>% # allows you to compute on a data frame a row-at-a-time
      dplyr::mutate(edge = paste0(sort(c(from, to)), collapse = "#")) %>%
      dplyr::group_by(edge) %>%
      dplyr::summarize(nb_strain_with_edge = n(),
                       #similarity = paste0(similarity, collapse = ",")) %>%
                       similarity = mean(similarity)) %>%
      ungroup() %>%
      tidyr::separate(edge, sep = "#", into = c("to", "from")) %>%
      dplyr::select(from, to, nb_strain_with_edge, similarity)

    all_node_trans_GSEA_KEGG <- rbind(emap_to_network_df(emap_Trans.5448.GSEA.KEGG.sig, what = "vertices"),
                                     emap_to_network_df(emap_Trans.HKU419.GSEA.KEGG.sig, what = "vertices"),
                                     emap_to_network_df(emap_Trans.PS003.GSEA.KEGG.sig, what = "vertices"),
                                     emap_to_network_df(emap_Trans.PS006.GSEA.KEGG.sig, what = "vertices"),
                                     emap_to_network_df(emap_Trans.SP444.GSEA.KEGG.sig, what = "vertices")) %>%
      dplyr::group_by(name) %>%
      dplyr::summarize(nb_strain_with_node = n(),
                       #gene_set_size = paste0(gene_set_size, collapse = ","),
                       gene_set_size = mean(gene_set_size),
                       p.adjust =paste0(p.adjust, collapse = ",")) %>%
      ungroup()

    # Check output
    head(all_edge_trans_GSEA_KEGG, 5)

    ## # A tibble: 5 × 4
    ##   from                            to              nb_strain_with_edge similarity
    ##   <chr>                           <chr>                         <int>      <dbl>
    ## 1 Fructose and mannose metabolism Amino sugar an…                   2      0.278
    ## 2 Galactose metabolism            Amino sugar an…                   1      0.261
    ## 3 Glycolysis / Gluconeogenesis    Amino sugar an…                   1      0.2  
    ## 4 Phosphotransferase system (PTS) Amino sugar an…                   2      0.288
    ## 5 Starch and sucrose metabolism   Amino sugar an…                   2      0.211

    head(all_node_trans_GSEA_KEGG, 5)

    ## # A tibble: 5 × 4
    ##   name                                nb_strain_with_node gene_set_size p.adjust
    ##   <chr>                                             <int>         <dbl> <chr>   
    ## 1 Amino sugar and nucleotide sugar m…                   2          18.5 0.00750…
    ## 2 Aminoacyl-tRNA biosynthesis                           1          15   0.04721…
    ## 3 Ascorbate and aldarate metabolism                     5          17.8 0.00980…
    ## 4 Biosynthesis of amino acids                           5          32.8 0.02125…
    ## 5 Biosynthesis of antibiotics                           2          52.5 0.00968…

    # 3. Save data
    write.table(all_edge_trans_GSEA_KEGG, "/Users/wmujchariyak/Desktop/Strep_all_edge_trans_GSEA_KEGG.tsv", sep = "\t", row.names = F)
    write.table(all_node_trans_GSEA_KEGG, "/Users/wmujchariyak/Desktop/Strep_all_node_trans_GSEA_KEGG.tsv", sep = "\t", row.names = F)
