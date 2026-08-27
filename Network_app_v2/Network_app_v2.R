#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)
library(visNetwork)  
library(svglite)     

# UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                titlePanel("Network of Enriched Pathways"),
                sidebarLayout(
                  sidebarPanel(
                    
                    # 1. Upload Data 
                    h5("1. Upload Data"),
                    fileInput("nodeFile", "Node Table (.tsv)", accept = ".tsv"),
                    fileInput("edgeFile", "Edge Table (.tsv)", accept = ".tsv"),
                    
                    # 2. Network Filter 
                    hr(),
                    h5("2. Network Filter"),
                    numericInput("numConditions", "Number of Conditions", 5, min=1, max=100, step=1),
                    radioButtons("filterType", "Filter Mode", 
                                 choices = c("Exact" = "exact", "Less Than or Equal" = "lte"),
                                 selected = "lte"),
                    actionButton("plotButton", "Plot Network"),
                    
                    # 3. Visualisation Style 
                    hr(),
                    h5("3. Visualisation Style"),
                    radioButtons("colorChoice", "Color Palette",
                                 choices = c("Preset" = "preset", "Custom Hex" = "custom")),
                    conditionalPanel(
                      condition = "input.colorChoice == 'preset'",
                      selectInput("paletteName", "Choose Preset Palette",
                                  choices = c("Blues", "Reds", "Viridis", "Plasma"))
                    ),
                    conditionalPanel(
                      condition = "input.colorChoice == 'custom'",
                      textInput("hexColors", "Enter Hex Codes (comma-separated):",
                                "#4575B4,#91BFDB,#FEE090,#FC8D59,#D73027")
                    ),
                    radioButtons("legendLabel", "Figure Legend Label",
                                 choices = c("Strains"    = "Strains",
                                             "Conditions" = "Conditions",
                                             "Omics"      = "Omics"),
                                 selected = "Conditions",
                                 inline   = TRUE),
                    
                    # 4. Interactive Plot 
                    hr(),
                    h5("4. Interactive Plot"),
                    selectInput("layoutType", "Layout",
                                choices = c(
                                  "Force-directed (default)" = "layout_with_fr",
                                  "Kamada-Kawai"            = "layout_with_kk",
                                  "Circle"                  = "layout_in_circle",
                                  "Hierarchical"            = "layout_as_tree"
                                ),
                                selected = "layout_with_kk"),
                    checkboxInput("showNodeLabels", "Show Node Labels", value = FALSE),
                    textInput("excludeNodes",
                              "Exclude Nodes (comma-separated names):",
                              placeholder = "e.g. NodeA, NodeB"),
                    tags$label("Node Grouping"),
                    checkboxInput("enableGrouping", "Enable node grouping", value = FALSE),
                    sliderInput("groupThreshold",
                                "Grouping Similarity Threshold (0–1):",
                                min = 0, max = 1, value = 0.8, step = 0.01),
                    helpText("Nodes connected by an edge with similarity ABOVE this threshold
                              will be merged into a single node. Names are joined with '/'."),
                    
                    
                    radioButtons("edgeThicknessVar", "Edge Thickness",
                                 choices = c("Similarity"               = "similarity",
                                             "Number of Conditions"  = "nb_condition_with_edge"),
                                 selected = "similarity",
                                 inline   = TRUE),
                    
                    # 5. Static Plot 
                    hr(),
                    h5("5. Static Plot"),
                    radioButtons("staticLayout", "Layout",
                                 choices = c("FR — clean overview"          = "fr",
                                             "KK — similarity-based clusters" = "kk"),
                                 selected = "kk"),
                    helpText(
                      tags$b("FR (Force-directed):"),
                      "Arranges nodes to minimise visual clutter.",
                      "Use for a clean overview of overall network structure.",
                      tags$br(),
                      tags$b("KK (Kamada-Kawai):"),
                      "Uses similarity scores to set distances. Similar pathways appear closer together."
                    ),
                    # 6. Node Ranking 
                    hr(),
                    h5("6. Node Ranking"),
                    radioButtons("rankColumn", "Rank Table By",
                                 choices = c("Interactions" = "Interactions",
                                             "Conditions"   = "Conditions"),
                                 selected = "Interactions",
                                 inline   = TRUE),
                    radioButtons("rankOrder", "Sort order:",
                                 choices = c("Highest first" = "desc",
                                             "Lowest first"  = "asc"),
                                 selected = "desc",
                                 inline   = TRUE),
                    numericInput("rankTop", "Show top N nodes in ranking:",
                                 value = 10, min = 1, max = 100, step = 1)
                    
                  ),
                  
                  mainPanel(
                    # Present Interactive and Static view
                    tabsetPanel(
                      tabPanel("Interactive",
                               visNetwork::visNetworkOutput("visNetworkPlot", height = "600px"),
                               helpText("Scroll to zoom, drag background to pan, drag nodes to rearrange."),
                               # Ranking table shown below the interactive plot
                               hr(),
                               h4("Node Interaction Ranking"),
                               tableOutput("nodeRankTable")
                               
                      ),
                      tabPanel("Static",
                               plotOutput("networkPlot"),
                               # Download buttons for the static figure
                               hr(),
                               h5("Download Static Figure:"),
                               downloadButton("download_png", "Download PNG"),
                               downloadButton("download_svg", "Download SVG")
                               
                      )
                    )
                    
                  )
                )
)

# Merge nodes whose connecting edge similarity >= threshold 
# Returns list(nodes, edges) with merged nodes named "A/B/C"
merge_similar_nodes <- function(nodes, edges, threshold) {
  
  # Find edges above threshold
  high_sim <- edges %>% filter(similarity >= threshold)
  
  if (nrow(high_sim) == 0) return(list(nodes = nodes, edges = edges))
  
  # Build a small graph from only the high-similarity edges to find connected components
  all_node_names <- unique(nodes$name)
  g_sim <- igraph::graph_from_data_frame(
    d        = high_sim[, c("from", "to")],
    vertices = data.frame(name = all_node_names),
    directed = FALSE
  )
  
  # Each connected component becomes one merged node
  comps      <- igraph::components(g_sim)
  membership <- data.frame(
    name       = names(comps$membership),
    cluster_id = as.integer(comps$membership),
    stringsAsFactors = FALSE
  )
  
  # Only merge components that have > 1 node
  comp_sizes   <- comps$csize
  solo_comps   <- which(comp_sizes == 1)
  merged_comps <- which(comp_sizes > 1)
  
  # Build the new merged node table
  new_nodes <- membership %>%
    left_join(nodes, by = "name") %>%
    group_by(cluster_id) %>%
    summarise(
      # merged name: sort alphabetically then join with "/"
      merged_name            = paste(sort(name), collapse = "/"),
      nb_condition_with_node = max(nb_condition_with_node),
      # union of all condition identifiers from merged nodes
      condition_info         = paste(unique(unlist(strsplit(condition_info, ","))), collapse = ","),
      # keep other numeric columns by taking the mean
      gene_set_size          = if ("gene_set_size" %in% names(.)) mean(gene_set_size, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    rename(name = merged_name) %>%
    # flag so we know which nodes are merged (for coloring if needed)
    mutate(is_merged = cluster_id %in% merged_comps) %>%
    select(-cluster_id)
  
  # Map original node names → new merged names for edge rewiring
  name_map <- membership %>%
    left_join(
      membership %>%
        left_join(nodes, by = "name") %>%
        group_by(cluster_id) %>%
        summarise(merged_name = paste(sort(name), collapse = "/"), .groups = "drop"),
      by = "cluster_id"
    ) %>%
    select(original = name, merged = merged_name)
  
  # Rewire edges: replace from/to with merged names, then drop internal edges
  new_edges <- edges %>%
    left_join(name_map, by = c("from" = "original")) %>%
    rename(from_merged = merged) %>%
    left_join(name_map, by = c("to" = "original")) %>%
    rename(to_merged = merged) %>%
    mutate(
      from = from_merged,
      to   = to_merged
    ) %>%
    select(-from_merged, -to_merged) %>%
    # Remove edges that are now within the same merged node (self-loops)
    filter(from != to) %>%
    # If multiple original edges now point to same merged pair, keep the one
    # with the highest nb_condition_with_edge (most representative)
    group_by(from, to) %>%
    slice_max(nb_condition_with_edge, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  list(nodes = new_nodes, edges = new_edges)
}

# Server
server <- function(input, output) {
  # Reactive to read TSV files
  nodeData <- reactive({
    req(input$nodeFile)
    read.delim(input$nodeFile$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  })
  edgeData <- reactive({
    req(input$edgeFile)
    read.delim(input$edgeFile$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  })
  # Reactive to create igraph object
  graphData <- eventReactive(input$plotButton, {
    nodes <- nodeData()
    edges <- edgeData()
    num <- input$numConditions
    
    # Filter nodes based on filter type
    if(input$filterType == "exact") {
      nodes_filtered <- nodes %>% filter(nb_condition_with_node == num)
    } else {
      nodes_filtered <- nodes %>% filter(nb_condition_with_node <= num)
    }
    
    # Filter edges that connect only to these nodes
    edge_nodes <- nodes_filtered$name
    edges_filtered <- edges %>%
      filter(from %in% edge_nodes & to %in% edge_nodes)
    
    ig <- igraph::graph_from_data_frame(d=edges_filtered, vertices=nodes_filtered, directed=FALSE)
    return(list(graph=ig, nodes=nodes_filtered, edges=edges_filtered))
    
  })
  
  # Reactive that optionally merges nodes based on similarity threshold
  groupedData <- reactive({
    req(graphData())
    nodes <- graphData()$nodes
    edges <- graphData()$edges
    
    if (isTRUE(input$enableGrouping)) {
      merge_similar_nodes(nodes, edges, input$groupThreshold)
    } else {
      list(nodes = nodes, edges = edges)
    }
  })
  
  # Interactive visNetwork plot
  output$visNetworkPlot <- visNetwork::renderVisNetwork({
    req(groupedData())
    nodes <- groupedData()$nodes
    edges <- groupedData()$edges
    
    
    # Exclude nodes by name
    excluded <- trimws(strsplit(input$excludeNodes, ",")[[1]])
    excluded <- excluded[excluded != ""]
    if (length(excluded) > 0) {
      nodes <- nodes %>% filter(!name %in% excluded)
      edges <- edges %>% filter(!from %in% excluded & !to %in% excluded)
    }
    # Build color palette 
    if (input$colorChoice == "preset") {
      palette_choice <- switch(input$paletteName,
                               "Blues"   = RColorBrewer::brewer.pal(9, "Blues"),
                               "Reds"    = RColorBrewer::brewer.pal(9, "Reds"),
                               "Viridis" = viridis::viridis(9),
                               "Plasma"  = viridis::plasma(9))
      color_palette <- colorRampPalette(palette_choice)(input$numConditions)
    } else {
      color_palette <- trimws(strsplit(input$hexColors, ",")[[1]])
    }
    
    # Prepare nodes for visNetwork (needs columns: id, label, color, value)
    vis_nodes <- nodes %>%
      mutate(
        id    = name,
        label = if (isTRUE(input$showNodeLabels)) name else "",  # toggle labels
        value = nb_condition_with_node,                   # controls node size
        color = color_palette[nb_condition_with_node],    # color by strain count
        title = paste0(
          "<b>", name, "</b><br>",
          input$legendLabel, ": ", nb_condition_with_node, "<br>",
          "<i>", gsub(",", ", ", condition_info), "</i>"
        )
      )
    
    # Edge thickness: use selected variable (similarity or nb_condition_with_edge)
    edge_var  <- if (input$edgeThicknessVar == "similarity") edges$similarity
    else edges$nb_condition_with_edge
    ev_range  <- range(edge_var, na.rm = TRUE)
    norm_ev   <- (edge_var - ev_range[1]) / max(ev_range[2] - ev_range[1], 1e-6)
    
    vis_edges <- edges %>%
      mutate(
        width = 1 + norm_ev * 7,                               # thickness encodes selected variable
        color = "#80808066",                                    # uniform grey
        title = paste0("Similarity: ",               round(similarity, 3),
                       "<br>Number of Conditions (edge): ", nb_condition_with_edge)  # hover tooltip
      )
    
    # Build the interactive plot
    visNetwork::visNetwork(vis_nodes, vis_edges) %>%
      visNetwork::visIgraphLayout(layout = input$layoutType) %>%
      visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visNetwork::visInteraction(navigationButtons = TRUE)
  })
  
  # Render node ranking table
  output$nodeRankTable <- renderTable({
    req(groupedData())
    nodes <- groupedData()$nodes
    edges <- groupedData()$edges
    
    # Count number of edges (interactions) per node
    edge_counts <- data.frame(name = c(edges$from, edges$to)) %>%
      group_by(name) %>%
      summarise(interactions = n(), .groups = "drop")
    
    # Join with node info and sort
    ranked <- nodes %>%
      left_join(edge_counts, by = "name") %>%
      mutate(interactions = replace_na(interactions, 0)) %>%
      select(Node = name, Conditions = nb_condition_with_node, Interactions = interactions)
    
    sort_col <- input$rankColumn   # "Interactions" or "Conditions"
    if (input$rankOrder == "desc") {
      ranked <- ranked %>% arrange(desc(.data[[sort_col]]))
    } else {
      ranked <- ranked %>% arrange(.data[[sort_col]])
    }
    
    head(ranked, input$rankTop)
  })
  
  output$networkPlot <- renderPlot({
    req(groupedData())
    nodes <- groupedData()$nodes
    edges <- groupedData()$edges
    
    # Exclude nodes by name (same as interactive plot)
    excluded <- trimws(strsplit(input$excludeNodes, ",")[[1]])
    excluded <- excluded[excluded != ""]
    if (length(excluded) > 0) {
      nodes <- nodes %>% filter(!name %in% excluded)
      edges <- edges %>% filter(!from %in% excluded & !to %in% excluded)
    }
    
    
    # Rebuild igraph from (optionally grouped) nodes and edges
    ig <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
    
    # 2. Add labels to the nodes 
    set.seed(1234) # for reproducibility
    tg <- tidygraph::as_tbl_graph(ig) %>% 
      tidygraph::activate(nodes) %>% 
      dplyr::mutate(label=name)
    
    # 3. color choice 
    if (input$colorChoice == "preset") {
      palette_choice <- switch(input$paletteName,
                               "Blues" = RColorBrewer::brewer.pal(9, "Blues"),
                               "Reds" = RColorBrewer::brewer.pal(9, "Reds"),
                               "Viridis" = viridis::viridis(9),
                               "Plasma" = viridis::plasma(9))
      color_palette <- colorRampPalette(palette_choice)(input$numConditions)
    } else {
      # Custom hex codes -> split by comma
      color_palette <- strsplit(input$hexColors, ",")[[1]]
      color_palette <- trimws(color_palette)  
    }
    
    
    
    # 4. Plot the network
    # Layout: FR (clean, unweighted) or KK (similarity-weighted, biological clusters)
    if (input$staticLayout == "kk") {
      edge_weights <- pmax(1 - edges$similarity, 1e-3)  # prevent zero weights crashing KK
      graph_layout <- ggraph::create_layout(tg, layout = "kk", weights = edge_weights)
    } else {
      graph_layout <- ggraph::create_layout(tg, layout = "fr")
    }
    

    # Determine which variable drives edge thickness and set legend label
    edge_width_var  <- if (input$edgeThicknessVar == "similarity") "similarity"
    else "nb_condition_with_edge"
    edge_width_label <- if (input$edgeThicknessVar == "similarity") "Similarity (edge)"
    else paste0("Number of ", input$legendLabel, " (edge)")
    
    
    plot_tg <- ggraph::ggraph(graph_layout) +
      
     
      # Edges: thickness encodes selected variable (similarity or nb_condition_with_edge)
      geom_edge_arc(color = "gray50", lineend = "round", strength = 0.1,
                    aes(edge_width = .data[[edge_width_var]], alpha = similarity)) +
   
      
      # Nodes: size + fill color both encode nb_strain (conserved = big + dark)
      ggraph::geom_node_point(
        aes(size  = nb_condition_with_node,
            color = nb_condition_with_node)
      ) +
      
      # Node labels
      geom_node_text(
        aes(label = name),
        repel         = TRUE,
        point.padding = unit(0.2, "lines"),
        color        = "gray10",
        size          = 3
      ) +
      
      # Node color scale → strain count
      scale_color_gradientn(
        colors = color_palette,
        limits = c(1, input$numConditions),
        breaks = 1:input$numConditions,
        name   = paste0("Number of ", input$legendLabel, "\n(node)")
      ) +
      
      # Node size scale → strain count
      scale_size_continuous(
        range  = c(2, input$numConditions * 1.5),
        limits = c(1, input$numConditions),
        breaks = 1:input$numConditions,
        name   = paste0("Number of ", input$legendLabel, "\n(node)")
      ) +
      
      
      # Edge width scale: label updates dynamically based on selected variable
      scale_edge_width(range = c(0.3, 3), name = edge_width_label,
                       guide = guide_legend(title.position = "top")) +

      scale_edge_alpha(range = c(0.2, 0.9), guide = "none") +
      
      theme_graph(background = "white") +
      theme(
        legend.position  = "bottom",
        legend.title     = element_text(size = 9, face = "bold"),
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.8, "lines"),
        legend.spacing.x = unit(0.6, "cm")
      ) +
      guides(
        color = guide_colorbar(
          title          = paste0("Number of ", input$legendLabel, " (node)"),
          title.position = "top",
          barwidth       = 8,
          barheight      = 0.5,
          order          = 1        # Node Color first
        ),
        size = guide_legend(
          title          = paste0("Number of ", input$legendLabel, " (node)"),
          title.position = "top",
          order          = 2        # Node Size second
        ),
        
        # Legend title for edge thickness updates dynamically
        edge_width = guide_legend(
          title          = edge_width_label,
          title.position = "top",
          order          = 3        # Edge Thickness third
        )
        
      )
    
    plot_tg
    
  })
  
  # Rebuilds the static plot — reused by both download handlers
  build_static_plot <- function() {
    req(groupedData())
    nodes <- groupedData()$nodes
    edges <- groupedData()$edges
    
    excluded <- trimws(strsplit(input$excludeNodes, ",")[[1]])
    excluded  <- excluded[excluded != ""]
    if (length(excluded) > 0) {
      nodes <- nodes %>% filter(!name %in% excluded)
      edges <- edges %>% filter(!from %in% excluded & !to %in% excluded)
    }
    
    ig <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
    
    set.seed(1234)
    tg <- tidygraph::as_tbl_graph(ig) %>%
      tidygraph::activate(nodes) %>%
      dplyr::mutate(label = name)
    
    if (input$colorChoice == "preset") {
      palette_choice <- switch(input$paletteName,
                               "Blues"   = RColorBrewer::brewer.pal(9, "Blues"),
                               "Reds"    = RColorBrewer::brewer.pal(9, "Reds"),
                               "Viridis" = viridis::viridis(9),
                               "Plasma"  = viridis::plasma(9))
      color_palette <- colorRampPalette(palette_choice)(input$numConditions)
    } else {
      color_palette <- trimws(strsplit(input$hexColors, ",")[[1]])
    }
    
    if (input$staticLayout == "kk") {
      edge_weights <- pmax(1 - edges$similarity, 1e-3)
      graph_layout <- ggraph::create_layout(tg, layout = "kk", weights = edge_weights)
    } else {
      graph_layout <- ggraph::create_layout(tg, layout = "fr")
    }
    

    # Determine which variable drives edge thickness and set legend label
    edge_width_var   <- if (input$edgeThicknessVar == "similarity") "similarity"
    else "nb_condition_with_edge"
    edge_width_label <- if (input$edgeThicknessVar == "similarity") "Similarity (edge)"
    else paste0("Number of ", input$legendLabel, " (edge)")

    
    ggraph::ggraph(graph_layout) +

      # Edges: thickness encodes selected variable
      geom_edge_arc(color = "gray50", lineend = "round", strength = 0.1,
                    aes(edge_width = .data[[edge_width_var]], alpha = similarity)) +

      ggraph::geom_node_point(
        aes(size  = nb_condition_with_node,
            color = nb_condition_with_node)
      ) +
      geom_node_text(
        aes(label = name),
        repel         = TRUE,
        point.padding = unit(0.2, "lines"),
        color        = "gray10",
        size          = 3
      ) +
      scale_color_gradientn(
        colors = color_palette,
        limits = c(1, input$numConditions),
        breaks = 1:input$numConditions,
        name   = paste0("Number of ", input$legendLabel, "\n(node)")
      ) +
      scale_size_continuous(
        range  = c(2, input$numConditions * 1.5),
        limits = c(1, input$numConditions),
        breaks = 1:input$numConditions,
        name   = paste0("Number of ", input$legendLabel, "\n(node)")
      ) +

      # Edge width scale: label updates dynamically based on selected variable
      scale_edge_width(range = c(0.3, 3), name = edge_width_label,
                       guide = guide_legend(title.position = "top")) +

      scale_edge_alpha(range = c(0.2, 0.9), guide = "none") +
      theme_graph(background = "white") +
      theme(
        legend.position  = "bottom",
        legend.title     = element_text(size = 9, face = "bold"),
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.8, "lines"),
        legend.spacing.x = unit(0.6, "cm")
      ) +
      guides(
        color = guide_colorbar(
          title          = paste0("Number of ", input$legendLabel, " (node)"),
          title.position = "top",
          barwidth       = 8,
          barheight      = 0.5
        ),
        size = guide_legend(
          title          = paste0("Number of ", input$legendLabel, " (node)"),
          title.position = "top"
        )
      )
  }
  
  # Download PNG
  output$download_png <- downloadHandler(
    filename = function() paste0("network_", input$staticLayout, "_", Sys.Date(), ".png"),
    content  = function(file) {
      ggplot2::ggsave(file, plot = build_static_plot(),
                      width = 14, height = 9, dpi = 300, device = "png")
    }
  )
  
  # Download SVG
  output$download_svg <- downloadHandler(
    filename = function() paste0("network_", input$staticLayout, "_", Sys.Date(), ".svg"),
    content  = function(file) {
      ggplot2::ggsave(file, plot = build_static_plot(),
                      width = 14, height = 9, device = svglite::svglite)
    }
  )
  
}
# Run the app
shinyApp(ui = ui, server = server)