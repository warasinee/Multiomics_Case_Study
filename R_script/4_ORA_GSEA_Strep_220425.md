ORA & GSEA
================

- [Load R packages](#load-r-packages)
- [Import all omics](#import-all-omics)
- [Set parameters](#set-parameters)
- [Example for only Transcriptomics](#example-for-only-transcriptomics)
- [For loop - ORA & GSEA analysis (All omics and all
  strains)](#for-loop---ora--gsea-analysis-all-omics-and-all-strains)
- [Save](#save)

## Load R packages

    library(clusterProfiler)

    ## 

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

    library(tidyverse)

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter()   masks clusterProfiler::filter(), stats::filter()
    ## ✖ dplyr::lag()      masks stats::lag()
    ## ✖ purrr::simplify() masks clusterProfiler::simplify()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

    library(conflicted)
    conflict_prefer("select", "dplyr")

    ## [conflicted] Will prefer dplyr::select over any other package.

    conflict_prefer("filter", "dplyr")

    ## [conflicted] Will prefer dplyr::filter over any other package.

    set.seed(123) # for reproducibility

## Import all omics

    # Import all entities fold change results (all omics combined)
    all_DE.df <- readRDS("/Users/wmujchariyak/Desktop/All_entity_fold_change_long_format.rds")

    # Import all functional annotations 
    all.annot <- readRDS("/Users/wmujchariyak/Desktop/All_funct_annot_clean_noduplicate.rda")

## Set parameters

    # For ORA only 
    my_FDR <- 0.05 # FDR threshold (filter signif. DE genes)
    my_logFC <- 1  # logFC threshold (filter signif. DE genes)

    # For ORA and GSEA
    annot_origin <- unique(all.annot$fterm_origin) %>% as.character() # Use all annotations
    is_core_accessory <- "yes"

## Example for only Transcriptomics

### ORA & GSEA

    ## Prepare Transcriptome data / K. pneumoniae
    my_species <- "Streptococcus pyogenes"
    my_omics <- "Transcriptomics RG"

    all_DE.res <- all_DE.df %>%
      filter(species == my_species) %>%
      filter(experiment == my_omics)

    if (grepl("Metabo", my_omics)){
      my.annot <- all.annot %>%
        filter(entity_type == "metabolite")
    } else {
      my.annot <- all.annot %>%
        filter(entity_type != "metabolite")
    }

    Trans.list.ORA <- list()
    Trans.list.ORA.summary <- list()

    Trans.list.GSEA <- list()
    Trans.list.GSEA.summary <- list()

    for (my_strain in unique(all_DE.res$strain)){
      # ORA analysis
      # prepare input data for ORA
      ORA.df <- all_DE.res %>%
        filter(is_species_core %in% is_core_accessory ) %>%
        filter(strain == my_strain) %>%
        filter(FDR < my_FDR) %>%
        filter(abs(logFC) > my_logFC) 
      
      glist.ORA <- ORA.df$entity_id
      
      # Get strain specific background gene list 
      bglist.ORA <- all_DE.res %>%
        filter(is_species_core %in% is_core_accessory) %>%
        filter(strain == my_strain) %>%
        distinct(entity_id, .keep_all = T) %>%
        .$entity_id
      
      term_to_gene.df <- my.annot %>%
        filter(fterm_origin %in% annot_origin) %>%
        filter(entity_id %in% bglist.ORA) %>%
        dplyr::select(fterm_id, entity_id) %>%
        distinct(fterm_id, entity_id) %>%
        ungroup()
      
      term_to_name.df <- my.annot %>%
        filter(fterm_origin %in% annot_origin) %>%
        filter(entity_id %in% bglist.ORA) %>%
        dplyr::select(fterm_id, fterm_name) %>%
        distinct(fterm_id, fterm_name) %>%
        ungroup()
      
      # Calculate mean/median logFC of term/path for plot
      ORA_logFC_summary <- ORA.df %>% 
        merge(., term_to_gene.df , by = "entity_id", all.x = T) %>%
        group_by(strain, fterm_id) %>%
        summarise(mean_logFC = mean(logFC), median_logFC = median(logFC)) %>%
        ungroup() %>%
        dplyr::select(-strain)
      
      ORA <- enricher(gene = glist.ORA,
                      TERM2GENE = term_to_gene.df, 
                      TERM2NAME = term_to_name.df,
                      universe = bglist.ORA,
                      qvalueCutoff = 1,
                      minGSSize = 5, 
                      maxGSSize = length(bglist.ORA)/2,
                      pvalueCutoff = 1,
                      pAdjustMethod = "fdr") 
      
      Trans.list.ORA[[my_strain]] <- ORA
      Trans.list.ORA.summary[[my_strain]]  <- ORA %>% 
        as.data.frame() %>% 
        mutate(strain = my_strain) %>%
        mutate(experiment = my_omics) %>%
        merge(., ORA_logFC_summary, by.x = "ID", by.y = "fterm_id") 
      
    ###############
      # GSEA analysis
      GSEA.df <- all_DE.res %>%
        filter(is_species_core %in% is_core_accessory ) %>%
        filter(strain == my_strain) %>%
        filter(experiment == my_omics) %>%
        distinct(entity_id, .keep_all = T) 
      
      glist.GSEA <- GSEA.df$logFC 
      names(glist.GSEA) <- GSEA.df$entity_id
      glist.GSEA <- sort(glist.GSEA, decreasing = T)
      
      GSEA <- GSEA(geneList = glist.GSEA, 
                   TERM2GENE = term_to_gene.df,
                   TERM2NAME = term_to_name.df,
                   minGSSize = 5,
                   maxGSSize = length(bglist.ORA)/2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr",
                   nPermSimple = 100000,  # Increase the number of permutations
                   eps = 0  # Allow more precise p-value estimation
                   )
      
      Trans.list.GSEA[[my_strain]] <- GSEA
      Trans.list.GSEA.summary[[my_strain]]  <- GSEA %>% 
        as.data.frame() %>% 
        mutate(strain = my_strain) %>%
        mutate(experiment = my_omics) %>%
        merge(., term_to_name.df, by.x = "ID", by.y = "fterm_id", all.x = T) 
    }

    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...

### Check ORA and GSEA Results

    head(Trans.list.ORA.summary[["5448"]], 5)

    ##           ID                      Description GeneRatio BgRatio    pvalue
    ## 1 GO:0000003                     reproduction     8/736 16/1450 0.6224925
    ## 2 GO:0000027 ribosomal large subunit assembly    12/736 22/1450 0.4437220
    ## 3 GO:0000028 ribosomal small subunit assembly     9/736 14/1450 0.2279206
    ## 4 GO:0000041   transition metal ion transport     6/736 16/1450 0.9068765
    ## 5 GO:0000049                     tRNA binding     9/736 25/1450 0.9553430
    ##    p.adjust    qvalue
    ## 1 0.9234277 0.8604017
    ## 2 0.8334554 0.7765702
    ## 3 0.6806156 0.6341620
    ## 4 0.9991873 0.9309905
    ## 5 0.9991873 0.9309905
    ##                                                                                                                                                                    geneID
    ## 1                                                         EW021_RS01715/EW021_RS02135/EW021_RS05600/EW021_RS06300/EW021_RS06380/EW021_RS06390/EW021_RS07640/EW021_RS08365
    ## 2 EW021_RS00385/EW021_RS00445/EW021_RS01585/EW021_RS02145/EW021_RS02150/EW021_RS03340/EW021_RS03410/EW021_RS03425/EW021_RS06420/EW021_RS07015/EW021_RS09115/EW021_RS09120
    ## 3                                           EW021_RS00380/EW021_RS00395/EW021_RS00435/EW021_RS01310/EW021_RS01585/EW021_RS03220/EW021_RS03475/EW021_RS03505/EW021_RS04910
    ## 4                                                                                     EW021_RS02115/EW021_RS02120/EW021_RS02125/EW021_RS03495/EW021_RS05980/EW021_RS08880
    ## 5                                           EW021_RS00430/EW021_RS02000/EW021_RS02150/EW021_RS02700/EW021_RS03330/EW021_RS03585/EW021_RS06450/EW021_RS08530/EW021_RS08550
    ##   Count strain         experiment mean_logFC median_logFC
    ## 1     8   5448 Transcriptomics RG -0.4393122    -1.121589
    ## 2    12   5448 Transcriptomics RG -1.8657007    -1.968675
    ## 3     9   5448 Transcriptomics RG -2.1933540    -2.395740
    ## 4     6   5448 Transcriptomics RG -1.9442218    -2.488583
    ## 5     9   5448 Transcriptomics RG -1.6087615    -1.583629

    head(Trans.list.GSEA.summary[["5448"]], 5)

    ##           ID                      Description setSize enrichmentScore
    ## 1 GO:0000003                     reproduction      16      -0.2558126
    ## 2 GO:0000027 ribosomal large subunit assembly      22      -0.6245913
    ## 3 GO:0000028 ribosomal small subunit assembly      14      -0.6339944
    ## 4 GO:0000041   transition metal ion transport      16      -0.5330212
    ## 5 GO:0000049                     tRNA binding      25      -0.5146465
    ##          NES       pvalue    p.adjust      qvalue rank
    ## 1 -0.7480318 0.8117691214 0.894599041 0.705741942  344
    ## 2 -1.9885974 0.0004347225 0.006859516 0.005411417  437
    ## 3 -1.7826521 0.0064415111 0.052657413 0.041541007  346
    ## 4 -1.5586287 0.0366229456 0.160488834 0.126608342  168
    ## 5 -1.6929170 0.0097196483 0.069418751 0.054763891  437
    ##                     leading_edge
    ## 1 tags=31%, list=24%, signal=24%
    ## 2 tags=77%, list=30%, signal=55%
    ## 3 tags=64%, list=24%, signal=49%
    ## 4 tags=31%, list=12%, signal=28%
    ## 5 tags=52%, list=30%, signal=37%
    ##                                                                                                                                                                                                                                 core_enrichment
    ## 1                                                                                                                                                                         EW021_RS08365/EW021_RS06380/EW021_RS06390/EW021_RS01715/EW021_RS07640
    ## 2 EW021_RS00370/EW021_RS04440/EW021_RS00365/EW021_RS00345/EW021_RS00350/EW021_RS02980/EW021_RS00385/EW021_RS03425/EW021_RS02150/EW021_RS02145/EW021_RS06420/EW021_RS07015/EW021_RS03340/EW021_RS01585/EW021_RS09120/EW021_RS09115/EW021_RS00445
    ## 3                                                                                                                 EW021_RS00380/EW021_RS00395/EW021_RS03505/EW021_RS01310/EW021_RS03475/EW021_RS01585/EW021_RS00435/EW021_RS03220/EW021_RS04910
    ## 4                                                                                                                                                                         EW021_RS02125/EW021_RS05980/EW021_RS02115/EW021_RS03495/EW021_RS02120
    ## 5                                                         EW021_RS00370/EW021_RS00345/EW021_RS00305/EW021_RS03390/EW021_RS08550/EW021_RS03585/EW021_RS02700/EW021_RS02000/EW021_RS08530/EW021_RS06450/EW021_RS02150/EW021_RS03330/EW021_RS00430
    ##   strain         experiment                       fterm_name
    ## 1   5448 Transcriptomics RG                     reproduction
    ## 2   5448 Transcriptomics RG ribosomal large subunit assembly
    ## 3   5448 Transcriptomics RG ribosomal small subunit assembly
    ## 4   5448 Transcriptomics RG   transition metal ion transport
    ## 5   5448 Transcriptomics RG                     tRNA binding

## For loop - ORA & GSEA analysis (All omics and all strains)

    # List of omics experiments
    omics_list <- c("Transcriptomics RG", "Proteomics MS1 DDA", "Metabolomics GC-MS")

    # Initialize lists to store results
    results_list <- list()

    # Loop through omics types
    for (my_omics in omics_list) {
      
      cat("\nProcessing Omics Experiment:", my_omics, "\n")
      
      # Prepare data for the current omics
      my_species <- "Streptococcus pyogenes"
      
      all_DE.res <- all_DE.df %>%
        filter(species == my_species) %>%
        filter(experiment == my_omics)
      
      # Filter annotation based on metabolite or non-metabolite
      if (grepl("Metabo", my_omics)) {
        my.annot <- all.annot %>%
          filter(entity_type == "metabolite")
      } else {
        my.annot <- all.annot %>%
          filter(entity_type != "metabolite")
      }
      
      # Initialize lists for ORA and GSEA results
      omics_list_ORA <- list()
      omics_list_ORA_summary <- list()
      
      omics_list_GSEA <- list()
      omics_list_GSEA_summary <- list()
      
      # Iterate over strains within the current omics experiment
      for (my_strain in unique(all_DE.res$strain)) {
        
        # ORA Analysis
        ORA.df <- all_DE.res %>%
          filter(is_species_core %in% is_core_accessory) %>%
          filter(strain == my_strain) %>%
          filter(FDR < my_FDR) %>%
          filter(abs(logFC) > my_logFC)
        
        glist.ORA <- ORA.df$entity_id
        
        # Background gene list
        bglist.ORA <- all_DE.res %>%
          filter(is_species_core %in% is_core_accessory) %>%
          filter(strain == my_strain) %>%
          distinct(entity_id, .keep_all = TRUE) %>%
          .$entity_id
        
        # Prepare term-to-gene and term-to-name mappings
        term_to_gene.df <- my.annot %>%
          filter(fterm_origin %in% annot_origin) %>%
          filter(entity_id %in% bglist.ORA) %>%
          dplyr::select(fterm_id, entity_id) %>%
          distinct(fterm_id, entity_id) %>%
          ungroup()
        
        term_to_name.df <- my.annot %>%
          filter(fterm_origin %in% annot_origin) %>%
          filter(entity_id %in% bglist.ORA) %>%
          dplyr::select(fterm_id, fterm_name) %>%
          distinct(fterm_id, fterm_name) %>%
          ungroup()
        
        # Calculate summary statistics for ORA
        ORA_logFC_summary <- ORA.df %>%
          merge(., term_to_gene.df, by = "entity_id", all.x = TRUE) %>%
          group_by(strain, fterm_id) %>%
          summarise(mean_logFC = mean(logFC), median_logFC = median(logFC)) %>%
          ungroup() %>%
          dplyr::select(-strain)
        
        # Perform ORA
        ORA <- enricher(
          gene = glist.ORA,
          TERM2GENE = term_to_gene.df,
          TERM2NAME = term_to_name.df,
          universe = bglist.ORA,
          qvalueCutoff = 1,
          minGSSize = 5,
          maxGSSize = length(bglist.ORA) / 2,
          pvalueCutoff = 1,
          pAdjustMethod = "fdr"
        )
        
        # Store ORA results
        omics_list_ORA[[my_strain]] <- ORA
        omics_list_ORA_summary[[my_strain]] <- ORA %>%
          as.data.frame() %>%
          mutate(strain = my_strain) %>%
          mutate(experiment = my_omics) %>%
          merge(., ORA_logFC_summary, by.x = "ID", by.y = "fterm_id")
        
        # GSEA Analysis
        GSEA.df <- all_DE.res %>%
          filter(is_species_core %in% is_core_accessory) %>%
          filter(strain == my_strain) %>%
          filter(experiment == my_omics) %>%
          distinct(entity_id, .keep_all = TRUE)
        
        glist.GSEA <- GSEA.df$logFC
        names(glist.GSEA) <- GSEA.df$entity_id
        glist.GSEA <- sort(glist.GSEA, decreasing = TRUE)
        
        GSEA <- GSEA(
          geneList = glist.GSEA,
          TERM2GENE = term_to_gene.df,
          TERM2NAME = term_to_name.df,
          minGSSize = 5,
          maxGSSize = length(bglist.ORA) / 2,
          pvalueCutoff = 1,
          pAdjustMethod = "fdr",
          nPermSimple = 100000,  # Increase the number of permutations
          eps = 0  # Allow more precise p-value estimation
        )
        
        # Store GSEA results
        omics_list_GSEA[[my_strain]] <- GSEA
        omics_list_GSEA_summary[[my_strain]] <- GSEA %>%
          as.data.frame() %>%
          mutate(strain = my_strain) %>%
          mutate(experiment = my_omics) %>%
          merge(., term_to_name.df, by.x = "ID", by.y = "fterm_id", all.x = TRUE)
      }
      
      # Save all results for this omics
      results_list[[my_omics]] <- list(
        ORA = omics_list_ORA,
        ORA_summary = omics_list_ORA_summary,
        GSEA = omics_list_GSEA,
        GSEA_summary = omics_list_GSEA_summary
      )
    }

    ## 
    ## Processing Omics Experiment: Transcriptomics RG

    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...

    ## 
    ## Processing Omics Experiment: Proteomics MS1 DDA

    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (19.69% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the `.groups` argument.preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (20.24% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the `.groups` argument.preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (17.39% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the `.groups` argument.preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (20.66% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the `.groups` argument.preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (20.13% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...
    ## done...

    ## 
    ## Processing Omics Experiment: Metabolomics GC-MS

    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...
    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the
    ## `.groups` argument.
    ## preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (0.78% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the `.groups` argument.preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (0.78% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...
    ## done...
    ## `summarise()` has grouped output by 'strain'. You can override using the `.groups` argument.preparing geneSet collections...
    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (1.55% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...
    ## done...

### Access results for a specific omics/strain

    # Access results for a specific omics
    # results_list[["Transcriptomics RG"]]

    # View ORA summary for 5448 (Transcriptome)
    head(results_list[["Transcriptomics RG"]]$ORA_summary[["5448"]], 5)

    ##           ID                      Description GeneRatio BgRatio    pvalue
    ## 1 GO:0000003                     reproduction     8/736 16/1450 0.6224925
    ## 2 GO:0000027 ribosomal large subunit assembly    12/736 22/1450 0.4437220
    ## 3 GO:0000028 ribosomal small subunit assembly     9/736 14/1450 0.2279206
    ## 4 GO:0000041   transition metal ion transport     6/736 16/1450 0.9068765
    ## 5 GO:0000049                     tRNA binding     9/736 25/1450 0.9553430
    ##    p.adjust    qvalue
    ## 1 0.9234277 0.8604017
    ## 2 0.8334554 0.7765702
    ## 3 0.6806156 0.6341620
    ## 4 0.9991873 0.9309905
    ## 5 0.9991873 0.9309905
    ##                                                                                                                                                                    geneID
    ## 1                                                         EW021_RS01715/EW021_RS02135/EW021_RS05600/EW021_RS06300/EW021_RS06380/EW021_RS06390/EW021_RS07640/EW021_RS08365
    ## 2 EW021_RS00385/EW021_RS00445/EW021_RS01585/EW021_RS02145/EW021_RS02150/EW021_RS03340/EW021_RS03410/EW021_RS03425/EW021_RS06420/EW021_RS07015/EW021_RS09115/EW021_RS09120
    ## 3                                           EW021_RS00380/EW021_RS00395/EW021_RS00435/EW021_RS01310/EW021_RS01585/EW021_RS03220/EW021_RS03475/EW021_RS03505/EW021_RS04910
    ## 4                                                                                     EW021_RS02115/EW021_RS02120/EW021_RS02125/EW021_RS03495/EW021_RS05980/EW021_RS08880
    ## 5                                           EW021_RS00430/EW021_RS02000/EW021_RS02150/EW021_RS02700/EW021_RS03330/EW021_RS03585/EW021_RS06450/EW021_RS08530/EW021_RS08550
    ##   Count strain         experiment mean_logFC median_logFC
    ## 1     8   5448 Transcriptomics RG -0.4393122    -1.121589
    ## 2    12   5448 Transcriptomics RG -1.8657007    -1.968675
    ## 3     9   5448 Transcriptomics RG -2.1933540    -2.395740
    ## 4     6   5448 Transcriptomics RG -1.9442218    -2.488583
    ## 5     9   5448 Transcriptomics RG -1.6087615    -1.583629

## Save

    #saveRDS(results_list, file = "/Users/wmujchariyak/Desktop/ORA_GSEA_list_Strep.rds")
