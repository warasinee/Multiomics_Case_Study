DIABLO results for network analysis
================

- [Load R packages](#load-r-packages)
- [Import all omics](#import-all-omics)
- [Data for MixOmics (DIABLO)](#data-for-mixomics-diablo)
- [Initial DIABLO model](#initial-diablo-model)
- [Tuning parameters](#tuning-parameters)
- [Final DIABLO model](#final-diablo-model)
- [Extract loading values](#extract-loading-values)
- [Import all omics (annotation)](#import-all-omics-annotation)
- [Set parameters](#set-parameters)
- [GSEA analysis](#gsea-analysis)
- [Network plot](#network-plot)
- [Network ananlysis to visualize interactions of enriched pathways](#network-ananlysis-to-visualize-interactions-of-enriched-pathways)

## Load R packages

    library(mixOmics)

    ## Loading required package: MASS

    ## Loading required package: lattice

    ## Loading required package: ggplot2

    ## 
    ## Loaded mixOmics 6.22.0
    ## Thank you for using mixOmics!
    ## Tutorials: http://mixomics.org
    ## Bookdown vignette: https://mixomicsteam.github.io/Bookdown
    ## Questions, issues: Follow the prompts at http://mixomics.org/contact-us
    ## Cite us:  citation('mixOmics')

    library(tidyverse)

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.4     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ✖ purrr::map()    masks mixOmics::map()
    ## ✖ dplyr::select() masks MASS::select()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

    library(igraph)

    ## 
    ## Attaching package: 'igraph'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     %--%, union
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union
    ## 
    ## The following objects are masked from 'package:purrr':
    ## 
    ##     compose, simplify
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     crossing
    ## 
    ## The following object is masked from 'package:tibble':
    ## 
    ##     as_data_frame
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     union

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
    ## 
    ## The following object is masked from 'package:igraph':
    ## 
    ##     simplify
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify
    ## 
    ## The following object is masked from 'package:lattice':
    ## 
    ##     dotplot
    ## 
    ## The following object is masked from 'package:MASS':
    ## 
    ##     select
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    library(RCy3)
    library(ggraph)

    set.seed(123) # for reproducibility

## Import all omics

    Strep_trans_data <- read.csv("/Users/wmujchariyak/Desktop/Strep_transcriptome_data.csv", header = TRUE) %>% column_to_rownames(var = "X")
    Strep_prot_data <- read.csv("/Users/wmujchariyak/Desktop/Strep_proteome_data.csv", header = TRUE) %>% column_to_rownames(var = "X")
    Strep_met_GC_data <- read.csv("/Users/wmujchariyak/Desktop/Strep_metabolome_GC_data.csv", header = TRUE) %>% column_to_rownames(var = "X")

## Data for MixOmics (DIABLO)

    data <- list(Transcript = Strep_trans_data, 
                 Protein = Strep_prot_data, 
                 Metabolite_GC = Strep_met_GC_data)
    lapply(data, dim) # check their dimensions (60 rows)

    ## $Transcript
    ## [1]   60 2239
    ## 
    ## $Protein
    ## [1]  60 993
    ## 
    ## $Metabolite_GC
    ## [1]  60 129

    Y <- factor(rep(c("RPMI", "SERUM"), 30)) # set the response variable as the Y df
    summary(Y)

    ##  RPMI SERUM 
    ##    30    30

## Initial DIABLO model

    design = matrix(0.1, ncol = length(data), nrow = length(data), 
                    dimnames = list(names(data), names(data)))
    diag(design) = 0 # set diagonal to 0s
    design

    ##               Transcript Protein Metabolite_GC
    ## Transcript           0.0     0.1           0.1
    ## Protein              0.1     0.0           0.1
    ## Metabolite_GC        0.1     0.1           0.0

## Tuning parameters

\*\* [LOO-CV
issue](https://mixomics-users.discourse.group/t/loo-cv-outputs-na-in-choice-ncomp/145/3) -
you can look at the perf plot and decide on the number of components
(Here, ncomp = 1).

    # Form basic DIABLO model with  an arbitrarily high number of components (ncomp = 5)
    plsda.diablo <- block.plsda(data, Y, ncomp = 5, design = design)

    ## Design matrix has changed to include Y; each block will be
    ##             linked to Y.

    set.seed(123) # for reproducibility
    # Tuning the number of components ("loo" = leave-one-out cross-validation, for small sample sizes)
    perf.plsda.diablo = perf(plsda.diablo, validation = 'loo', nrepeat = 1) 

    #perf.plsda.diablo$error.rate  # Lists the different types of error rates
    # Plot of the error rates based on weighted vote
    plot(perf.plsda.diablo) # plot output of tuning

![](7_DIABLO_Network_Strep_files/figure-markdown_strict/unnamed-chunk-5-1.png)

## Final DIABLO model

    # Here, we used block.plsda which is different from block.splsda (no variable selection)
    final.plsda.diablo <- block.plsda(data, Y, ncomp = 1, design = design)

    ## Design matrix has changed to include Y; each block will be
    ##             linked to Y.

## Extract loading values 

    # We used loading values for ranking importanct viables
    loadings_data <- final.plsda.diablo$loadings 

    loadings_data_trans <- loadings_data$Transcript
    loadings_data_prot <- loadings_data$Protein
    loadings_data_met <- loadings_data$Metabolite_GC

## Import all omics (annotation)

    # Import all entities fold change results (all omics combined)
    all_DE.df <- readRDS("/Users/wmujchariyak/Desktop/All_entity_fold_change_long_format.rds")
    all_DE.df.short <- all_DE.df %>%  dplyr::select(entity_id, core_entity_id) %>% unique()

    # Import all functional annotations 
    all.annot <- readRDS("/Users/wmujchariyak/Desktop/All_funct_annot_clean_noduplicate.rda")

## Set parameters

    # For ORA and GSEA
    annot_origin <- unique(all.annot$fterm_origin) %>% as.character() # Use all annotations

## GSEA analysis

    ## Prepare Transcriptome data
    my_species <- "Streptococcus pyogenes"
    my_omics <- "Transcriptomics RG"

    all_DE.res <- all_DE.df %>%
      filter(species == my_species) %>%
      filter(experiment == my_omics)

    # Annotation - add core_entity_id 
    if (grepl("Metabo", my_omics)){
      my.annot <- all.annot %>%
        filter(entity_type == "metabolite")
    } else {
      my.annot <- all.annot %>%
        filter(entity_type != "metabolite")
    }


    # Loading values 
    GSEA.df.trans <- loadings_data_trans %>% as.data.frame(.) %>% 
                     rownames_to_column(., var = "core_entity_id") %>% 
                     rename(loadings = comp1)

    # Get background gene list  (all variables from MixOmics)
    bglist.trans <- GSEA.df.trans$core_entity_id

    term_to_gene.df.trans <- my.annot %>%
        filter(fterm_origin %in% annot_origin) %>% 
        left_join(., all_DE.df.short, by = "entity_id") %>% 
        filter(core_entity_id %in% bglist.trans) %>%
        distinct(fterm_id, core_entity_id) %>%
        ungroup()

    ## Warning in left_join(., all_DE.df.short, by = "entity_id"): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 625 of `x` matches multiple rows in `y`.
    ## ℹ Row 63083 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

    term_to_name.df.trans <- my.annot %>%
        filter(fterm_origin %in% annot_origin) %>%
        left_join(., all_DE.df.short, by = "entity_id") %>% 
        filter(core_entity_id %in% bglist.trans) %>%
        dplyr::select(fterm_id, fterm_name) %>%
        distinct(fterm_id, fterm_name) %>%
        ungroup()

    ## Warning in left_join(., all_DE.df.short, by = "entity_id"): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 625 of `x` matches multiple rows in `y`.
    ## ℹ Row 63083 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

    #note: Warning: Detected an unexpected many-to-many relationship between `x` and `y`.

    ###########################################################
    ## GSEA analysis

    glist.GSEA.trans <- GSEA.df.trans$loadings
    names(glist.GSEA.trans) <- GSEA.df.trans$core_entity_id
    glist.GSEA.trans <- sort(glist.GSEA.trans, decreasing = T)

    GSEA.output.trans <- GSEA(geneList = glist.GSEA.trans, 
                   TERM2GENE = term_to_gene.df.trans,
                   TERM2NAME = term_to_name.df.trans,
                   minGSSize = 5,
                   maxGSSize = length(bglist.trans)/2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr",
                   nPermSimple = 100000,  # Increase the number of permutations
                   eps = 0  # Allow more precise p-value estimation
                   )

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

# Save

    #saveRDS(GSEA.output.trans, file = "/Users/wmujchariyak/Desktop/Strep.DIABLO.GSEA.output.trans.rds")

## Network plot 

    # Example - DIABLO.GSEA.KEGG.sig
    # 1. Prepare data
    Trans.DIABLO.GSEA.KEGG.sig <- clusterProfiler::filter(GSEA.output.trans, p.adjust < 0.05, grepl("^map", ID))

    # 2. Extract node and edge data from Trans.GSEA.KEGG 
    # 2.1 Calculate pairwise similarity (of the enriched terms using Jaccard’s similarity index)
          # Check added similarity matrix in the termsim slot of enrichment result 
          # For example, View(Trans.DIABLO.GSEA.KEGG.sig_2@termsim)
    Trans.DIABLO.GSEA.KEGG.sig_2 <- enrichplot::pairwise_termsim(Trans.DIABLO.GSEA.KEGG.sig, method = "JC")  

    # 2.2 Create Enrichment map 
    emap_Trans.DIABLO.GSEA.KEGG.sig <- enrichplot::emapplot(Trans.DIABLO.GSEA.KEGG.sig_2)

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

    # 2.4 Prepare edges df and node df of 1 strain (Please review 5_Node_Edge_data_Strep_220425.md for multiple strains)
    all_edge_DIABLO_trans_GSEA_KEGG <- emap_to_network_df(emap_Trans.DIABLO.GSEA.KEGG.sig) %>%
      dplyr::rowwise() %>% # allows you to compute on a data frame a row-at-a-time
      dplyr::mutate(edge = paste0(sort(c(from, to)), collapse = "#")) %>%
      dplyr::group_by(edge) %>%
      dplyr::summarize(nb_strain_with_edge = n(),
                       #similarity = paste0(similarity, collapse = ",")) %>%
                       similarity = mean(similarity)) %>%
      ungroup() %>%
      tidyr::separate(edge, sep = "#", into = c("to", "from")) %>%
      dplyr::select(from, to, nb_strain_with_edge, similarity)


    all_node_DIABLO_trans_GSEA_KEGG <- emap_to_network_df(emap_Trans.DIABLO.GSEA.KEGG.sig, what = "vertices")  %>%
      dplyr::group_by(name) %>%
      dplyr::summarize(nb_strain_with_node = n(),
                       #gene_set_size = paste0(gene_set_size, collapse = ","),
                       gene_set_size = mean(gene_set_size),
                       p.adjust =paste0(p.adjust, collapse = ",")) %>%
                       ungroup()

# save

    #write.table(all_edge_DIABLO_trans_GSEA_KEGG, "/Users/wmujchariyak/Desktop/Strep_all_edge_DIABLO_trans_GSEA_KEGG.tsv", sep = "\t", row.names = F)
    #write.table(all_node_DIABLO_trans_GSEA_KEGG, "/Users/wmujchariyak/Desktop/Strep_all_node_DIABLO_trans_GSEA_KEGG.tsv", sep = "\t", row.names = F)

## Network ananlysis to visualize interactions of enriched pathways

    # Example for only transcriptomics (GSEA - KEGG)

    # 1. Create igraph graphs from data frames
    # Check: class(ig_trans_GSEA_KEGG)
    ig_trans_GSEA_KEGG <- igraph::graph_from_data_frame(d=all_edge_DIABLO_trans_GSEA_KEGG, vertices=all_node_DIABLO_trans_GSEA_KEGG, directed = FALSE)

    # 2. Add labels to the nodes 
    tg_trans_GSEA_KEGG <- tidygraph::as_tbl_graph(ig_trans_GSEA_KEGG) %>% 
      tidygraph::activate(nodes) %>% 
      dplyr::mutate(label=name)

    # 3. Plot the network
    # set seed (so that the exact same network graph is created every time)
    set.seed(12345)

    # 3.1 Create a graph using the 'tg' data frame with the Fruchterman-Reingold layout
    plot_tg_trans_GSEA_KEGG_3 <- tg_trans_GSEA_KEGG %>%
      ggraph::ggraph(layout = "fr") +
      
      # Add arcs for edges with various aesthetics
      geom_edge_arc(colour = "gray50",
                    lineend = "round",
                    strength = .1,
                    aes(edge_width = all_edge_DIABLO_trans_GSEA_KEGG$nb_strain_with_edge,
                        alpha = all_edge_DIABLO_trans_GSEA_KEGG$nb_strain_with_edge)) +
      
      # Add points for nodes with size and color (for 1 strain)
      ggraph::geom_node_point(size = 5,
                              aes(color = all_node_DIABLO_trans_GSEA_KEGG$nb_strain_with_node)) +
      
      # Add edge of nodes on the next layer (shape = 21)
      ggraph::geom_node_point(size = 5, shape = 21) + 
      
      # Add text labels for nodes with various aesthetics
      geom_node_text(aes(label = name), 
                     repel = TRUE, 
                     point.padding = unit(0.2, "lines"), 
                     colour = "gray10") +
      
      # Change color palette with custom colors
      scale_color_gradientn(
        colors ="#D73027") +
      
      # Adjust edge width and alpha scales
      scale_edge_width(range = c(0.5, 2.5)) +
      scale_edge_alpha(range = c(0.2, 0.5)) +
      
      # Set graph background color to white
      theme_graph(background = "white") +
      
      # Adjust legend position to the top
      theme(legend.position = "bottom", 
            # suppress legend title
            legend.title = element_blank()) 

    # 3.2 View the plot 
    plot_tg_trans_GSEA_KEGG_3 # with fig.height=8, fig.width=20

![](7_DIABLO_Network_Strep_files/figure-markdown_strict/unnamed-chunk-15-1.png)
