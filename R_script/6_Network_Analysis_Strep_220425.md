Network Analysis - Visualize interactions of enriched pathways
================

- [Load R packages](#load-r-packages)
- [Import node and edge data](#import-node-and-edge-data)
- [Network ananlysis to visualize interactions of enriched
  pathways](#network-ananlysis-to-visualize-interactions-of-enriched-pathways)

## Load R packages

    library(igraph)

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

    library(ggraph)

    ## Loading required package: ggplot2

    library(tidyverse)

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.4     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%--%()      masks igraph::%--%()
    ## ✖ dplyr::as_data_frame() masks tibble::as_data_frame(), igraph::as_data_frame()
    ## ✖ purrr::compose()       masks igraph::compose()
    ## ✖ tidyr::crossing()      masks igraph::crossing()
    ## ✖ dplyr::filter()        masks stats::filter()
    ## ✖ dplyr::lag()           masks stats::lag()
    ## ✖ purrr::simplify()      masks igraph::simplify()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

    library(tidygraph)

    ## 
    ## Attaching package: 'tidygraph'
    ## 
    ## The following object is masked from 'package:igraph':
    ## 
    ##     groups
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    set.seed(12345) # for reproducibility

## Import node and edge data

    ## Transcriptome
    all_edge_df_trans_GSEA_KEGG <- read.table("/Users/wmujchariyak/Desktop/Strep_all_edge_trans_GSEA_KEGG.tsv", header = TRUE)
    all_node_df_trans_GSEA_KEGG <- read.table("/Users/wmujchariyak/Desktop/Strep_all_node_trans_GSEA_KEGG.tsv", header = TRUE)

## Network ananlysis to visualize interactions of enriched pathways

    # Example for only transcriptomics (GSEA - KEGG)

    # 1. Create igraph graphs from data frames
    # Check: class(ig_trans_GSEA_KEGG)
    ig_trans_GSEA_KEGG <- igraph::graph_from_data_frame(d=all_edge_df_trans_GSEA_KEGG, vertices=all_node_df_trans_GSEA_KEGG, directed = FALSE)

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
                    aes(edge_width = all_edge_df_trans_GSEA_KEGG$nb_strain_with_edge,
                        alpha = all_edge_df_trans_GSEA_KEGG$nb_strain_with_edge)) +
      
      # Add points for nodes with size based on nb_strain_with_node and color based on gene_set_size
      ggraph::geom_node_point(size = all_node_df_trans_GSEA_KEGG$nb_strain_with_node*1.5,
                              aes(color = all_node_df_trans_GSEA_KEGG$nb_strain_with_node)) +
      
      # Add edge of nodes on the next layer (shape = 21)
      ggraph::geom_node_point(size = all_node_df_trans_GSEA_KEGG$nb_strain_with_node*1.5, shape = 21) + 
      
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
    plot_tg_trans_GSEA_KEGG_3 # with fig.height=8, fig.width=20

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/6_Network_plot_1.png)
