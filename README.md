# Multi-omics: Case Study 

This is a case study using the multi-omics analysis pipeline (XXX). 

In this study, we applied multivariate data integration methods and network analysis to previously published multi-omics datasets [(Mu et al., 2023)](https://www.nature.com/articles/s41467-023-37200-w) generated from five *Streptococcus pyogenes* and two *Klebsiella pneumoniae* genotypes exposed to human serum. 

## Case study 1 - *Streptococcus pyogenes*
## 2. Multi-omics Integration by MixOmics  <img src="https://github.com/warasinee/Multiomics_Analyses_2024/blob/main/Image/MixOmics_Logo.png" width=5% height=5%>

We implemented the MixOmics data analytic pipelines including PCA and Multiblock (s)PLS-DA or DIABLO. 

### **[2.1. PCA](https://github.com/warasinee/Multiomics_Analyses_2024/blob/main/R_script/1_PCA_130125.md)** 
Initially, Principal Component Analysis (PCA), a dimensionality reduction and unsupervised machine learning method, was performed to assess the similarities of bacterial responses to distinct media conditions (Serum and RPMI).  

**Key steps:** 
```r
library(mixOmics)
set.seed(123) # for reproducibility

data <- list(Transcript = trans_data, 
             Protein = prot_data, 
             Metabolite_GC = met_GC_data, 
             Metabolite_LC = met_LC_data)
lapply(data, dim) # check their dimensions
Y <- factor(rep(c("RPMI", "SERUM"), 20)) # set the response variable as the Y df
summary(Y)

data_name <- c("Transcript", "Protein", "Metabolite_GC", "Metabolite_LC")
list.final.pca <- list()
plot.tune.pca <- list()
pca.cum.var <- list()

# PCA loop
for (i in seq_along(data_name)){
  data_pca <- data[[i]]
  # Choosing number of components
  tune.pca <- tune.pca(data_pca, ncomp = 10, scale = TRUE)
  plot.tune.pca[[i]] <- tune.pca
  # Outputs cumulative proportion of variance
  pca.cum.var[[i]] <- tune.pca$cum.var       
  # Final pca 
  final.pca <- pca(data_pca, ncomp = 3, scale = TRUE)
  list.final.pca[[i]] <- final.pca
}
# Note: Results from "tune.pca" suggest that using 3 components can cover almost 50% of variation in the datasets
# Visualization 
plotIndiv_pca_trans <- plotIndiv(list.final.pca[[1]], group = Y, ind.names = FALSE,
                                  legend = TRUE, ellipse = T, style = "ggplot2", size.title = rel(2.5),
                                  size.xlabel = rel(1.5), size.ylabel = rel(1.5),
                                  X.label = "PC1: 24%",
                                  Y.label = "PC2: 15%",  
                                  title = 'Transcriptome') 
```

### **[2.2. DIABLO](https://github.com/warasinee/Multiomics_Analyses_2024/blob/main/R_script/2_DIABLO_130125.md)** 
Here, we performed integrative supervised analysis using Data Integration Analysis for Biomarker discovery using Latent cOmponents (DIABLO), providing highly correlated variables from omics datasets discriminating *S. aureus* responses under two different conditions (Serum and RPMI).
This script contains all steps to identify signature molecules of *S. aureus* serum responses, including generating basic DIABLO model, tuning the number of components, and creating final DIABLO model.

**Key steps:** 
```r
library(mixOmics)
set.seed(123) # for reproducibility

# Data for MixOmics
data <- list(Transcript = trans_data, 
             Protein = prot_data, 
             Metabolite_GC = met_GC_data, 
             Metabolite_LC = met_LC_data)
lapply(data, dim) # check their dimensions
Y <- factor(rep(c("RPMI", "SERUM"), 20)) # set the response variable as the Y df
summary(Y)
#######################################
# Initial DIABLO model
design = matrix(0.5, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s
design
#######################################
# Tuning parameters - number of components and number of features 
# Form basic DIABLO model with  an arbitrarily high number of components (ncomp = 5)
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design)

## 1. Tuning the number of components ("loo" = leave-one-out cross-validation, for small sample sizes)
perf.diablo = perf(basic.diablo.model, validation = 'loo', nrepeat = 1) 
# Plot output of tuning to select ncomp
plot(perf.diablo) # select ncomp = 1

## 2. Tuning the number of features 
** This may take some time to run (~ 15-20 min)

# Set grid of values for each component to test
test.keepX = list (Transcript = seq(20,100,20), 
                   Protein = seq(20,100,20),
                   Metabolite_GC = seq(10,100,10),
                   Metabolite_GC = seq(10,100,10)) # Note: seq(10,100,10) = start from 10 and increase by 10 until 100

# Run the feature selection tuning (ncomp = 1 - from tuning the number of components)
tune.SA = tune.block.splsda(X = data, Y = Y, ncomp = 1, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 5, nrepeat = 10,
                              dist = "centroids.dist")

# Set the optimal values of features to retain
list.keepX = tune.SA$choice.keepX 
list.keepX 
#######################################
# Final DIABLO model
# Final model (use arbitrary number for # variables from our preliminary results in "Tuning the number of features")

list.keepX <-  list (Transcript =  20, 
                    Protein =  20,
                    Metabolite_GC = 10,
                    Metabolite_LC =  10) 

# Set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = 1, 
                                  keepX = list.keepX, design = design)

final.diablo.model$design # design matrix for the final model
#######################################
# AUC for the model (Example)
auc.diablo.trans.com <- auroc(final.diablo.model, roc.block = "Transcript", roc.comp = 1,
                        print = FALSE)
#######################################
# Visualization
## 1. Sample plots
plotDiablo(final.diablo.model, ncomp = 1)
## 2. Plot Loadings
plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median', size.name = 0.65)
## 3. Circos plot 
circosPlot(final.diablo.model, cutoff = 0.9, comp = 1, line = TRUE, 
           color.blocks = c('darkorchid', 'brown1', 'lightgreen',"orange"),
           color.cor = c("chocolate3","grey20"), size.labels = 1.2, size.variables = 0.5, size.legend = 1)
#######################################
# Extract data from final DIABLO model
## 1. Variables selected on component 1
trans_var <- selectVar(final.diablo.model, block = 'Transcript', comp = 1)
## 2. Correlation matrix from the circos plot
corMat <- circosPlot(final.diablo.model, cutoff = 0.7)
corMat2 <- as.data.frame(corMat)
#######################################
# Save data for further visualization in Cytoscape
library(igraph)
myNetwork <- network(final.diablo.model, blocks = c(1,2,3,4), cutoff = 0.9) 
write_graph(myNetwork$gR, file = "/Users/wmujchariyak/Desktop/myNetwork_conserved.gml", format = "gml")
```

### **[2.3. Cytoscape](https://github.com/warasinee/Multiomics_Analyses_2024/blob/main/R_script/3_Cytoscape_121224.R)** 
In this step, you need to connect R to Cytoscape software using [RCy3 package](https://cytoscape.org/RCy3/articles/Overview-of-RCy3.html).

**Key steps:** 
```r
# Import .gml file (output from MixOmics (DIABLO))
library(igraph)
my.network.conserved <- igraph::read_graph("Path_to_file/myNetwork_conserved.gml",format=c("gml"))

# From igraph to Cytoscape (Now you will see the network in Cytoscape)
library(RCy3)
createNetworkFromIgraph(my.network.conserved,"myIgraph.conserved")
```

## 3. Network and Enrichment Analyses
We performed Pathway Enrichment Analysis (PEA) including both an overrepresentation analysis (ORA) and a gene set enrichment analysis (GSEA), to reduce the complexity of data and discern the overrepresented biological pathways.
### **[3.1. ORA/GESA](https://github.com/warasinee/Multiomics_Analyses_2024/blob/main/R_script/4_ORA_GSEA_140125.md)**
**Key steps:** 
```r
library(clusterProfiler)
set.seed(123) # for reproducibility

## Prepare Transcriptome data / S.aureus
my_species <- "Staphylococcus aureus"
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
#######################################
## ORA & GSEA loop
 
Trans.list.ORA <- list()
Trans.list.ORA.summary <- list()
Trans.list.GSEA <- list()
Trans.list.GSEA.summary <- list()

for (my_strain in unique(all_DE.res$strain)){
  # ORA analysis
  # Prepare input data for ORA
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
  
#######################################
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
               eps = 0  # Allow more precise p-value estimation)
   Trans.list.GSEA[[my_strain]] <- GSEA
  Trans.list.GSEA.summary[[my_strain]]  <- GSEA %>% 
    as.data.frame() %>% 
    mutate(strain = my_strain) %>%
    mutate(experiment = my_omics) %>%
    merge(., term_to_name.df, by.x = "ID", by.y = "fterm_id", all.x = T) 
}
#######################################
# Check the results 
head(Trans.list.ORA.summary[["BPH2760"]], 5)
head(Trans.list.GSEA.summary[["BPH2760"]], 5)
```
### **[3.2. Node & Edge data](https://github.com/warasinee/Multiomics_Analyses_2024/blob/main/R_script/5_Node_Edge_data_140125.md)**
To prepare node and edge data for network analysis, we used pairewise_termsim function in the enrichplot package in R. This function calculates the pairwise similarity of the enriched terms using Jaccard’s similarity coefficient or the similarity of gene subsets sharing between pathways.

Then we wrote the function in R to convert enrichment map to dataframes which provide similarity scores, and number of strains shared individual enriched pathways representing node and edge, respectively, for the network visualization.     

**Key steps:**  
```r
# Prepare node and edge data for network anlysis
# 1. Prepare data - select only significant pathways (p.adjust < 0.05)
  ## For example, subseting of significant (KEGG) pathways from ORA results
  Trans.BPH2760.ORA.KEGG.sig <- clusterProfiler::filter(Trans.list.ORA$BPH2760, p.adjust < 0.05, grepl("^map", ID))

# 2. Extract node and edge data from ORA and GSEA results
  # 2.1 Calculate pairwise similarity
  ## For example, calculating pairwise similarity of significant (KEGG) pathways
  Trans.BPH2760.ORA.KEGG.sig_2 <- enrichplot::pairwise_termsim(Trans.BPH2760.ORA.KEGG.sig, method = "JC")
 # 2.2 Create enrichment map
  ## For example, creating enrichment map from the previous step
  emap_Trans.BPH2760.ORA.KEGG.sig <- enrichplot::emapplot(Trans.BPH2760.ORA.KEGG.sig_2)

  # 2.3 Function to converse enrichment map to dataframe(df)
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
 # 2.4 Concatenate node df and edge df of the 5 strains
all_edge_trans_ORA_KEGG <- rbind(emap_to_network_df(emap_Trans.BPH2760.ORA.KEGG.sig),
                                 emap_to_network_df(emap_Trans.BPH2819.ORA.KEGG.sig),
                                 emap_to_network_df(emap_Trans.BPH2900.ORA.KEGG.sig),
                                 emap_to_network_df(emap_Trans.BPH2947.ORA.KEGG.sig),
                                 emap_to_network_df(emap_Trans.BPH2986.ORA.KEGG.sig)) %>%
  dplyr::rowwise() %>% # allows you to compute on a data frame a row-at-a-time
  dplyr::mutate(edge = paste0(sort(c(from, to)), collapse = "#")) %>%
  dplyr::group_by(edge) %>%
  dplyr::summarize(nb_strain_with_edge = n(),
                   #similarity = paste0(similarity, collapse = ",")) %>%
                   similarity = mean(similarity)) %>%
  ungroup() %>%
  tidyr::separate(edge, sep = "#", into = c("to", "from")) %>%
  dplyr::select(from, to, nb_strain_with_edge, similarity)

all_node_trans_ORA_KEGG <- rbind(emap_to_network_df(emap_Trans.BPH2760.ORA.KEGG.sig, what = "vertices"),
                                 emap_to_network_df(emap_Trans.BPH2819.ORA.KEGG.sig, what = "vertices"),
                                 emap_to_network_df(emap_Trans.BPH2900.ORA.KEGG.sig, what = "vertices"),
                                 emap_to_network_df(emap_Trans.BPH2947.ORA.KEGG.sig, what = "vertices"),
                                 emap_to_network_df(emap_Trans.BPH2986.ORA.KEGG.sig, what = "vertices")) %>%
  dplyr::group_by(name) %>%
  dplyr::summarize(nb_strain_with_node = n(),
                   #gene_set_size = paste0(gene_set_size, collapse = ","),
                   gene_set_size = mean(gene_set_size),
                   p.adjust =paste0(p.adjust, collapse = ",")) %>% ungroup()

# Check output
head(all_edge_trans_ORA_KEGG, 5)
head(all_node_trans_ORA_KEGG, 5) 
```
### **[3.3. Network Visualization](https://github.com/warasinee/Multiomics_Analyses_2024/blob/main/R_script/6_Network_Analysis_140125.md)** 
Now we can visualize the network using node and edge data from previous step. 

**Key steps:**
```r
set.seed(12345) # for reproducibility
# 1. Create igraph graph from dataframe
ig_trans_ORA_KEGG <- igraph::graph_from_data_frame(d=all_edge_df_trans_ORA_KEGG, vertices=all_node_df_trans_ORA_KEGG, directed = FALSE)

# 2. Add lables to nodes
tg_trans_ORA_KEGG <- tidygraph::as_tbl_graph(ig_trans_ORA_KEGG) %>% 
  tidygraph::activate(nodes) %>% 
  dplyr::mutate(label=name)

# 3. Plot the network
  # 3.1 Create a graph using the 'tg' data frame with the Fruchterman-Reingold layout
plot_tg_trans_ORA_KEGG_3 <- tg_trans_ORA_KEGG %>%
  ggraph::ggraph(layout = "fr") +
  
  # Add arcs for edges with various aesthetics
  geom_edge_arc(colour = "gray50",
                lineend = "round",
                strength = .1,
                aes(edge_width = all_edge_df_trans_ORA_KEGG$nb_strain_with_edge,
                    alpha = all_edge_df_trans_ORA_KEGG$nb_strain_with_edge)) +
  
  # Add points for nodes with size based on nb_strain_with_node and color based on gene_set_size
  ggraph::geom_node_point(size = all_node_df_trans_ORA_KEGG$nb_strain_with_node*1.5,
                          aes(color = all_node_df_trans_ORA_KEGG$nb_strain_with_node)) +
  
  # Add edge of nodes on the next layer (shape = 21)
  ggraph::geom_node_point(size = all_node_df_trans_ORA_KEGG$nb_strain_with_node*1.5, shape = 21) + 
   # Add text labels for nodes with various aesthetics
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 point.padding = unit(0.2, "lines"), 
                 colour = "gray10") +
  
  # Change color palette with five custom colors
  scale_color_gradientn(
    colors = c("#4575B4", "#91BFDB", "#FEE090", "#FC8D59", "#D73027")) +
  
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
plot_tg_trans_ORA_KEGG_3 # with fig.height=8, fig.width=20
```





## Case study 2 - *Klebsiella pneumoniae*
