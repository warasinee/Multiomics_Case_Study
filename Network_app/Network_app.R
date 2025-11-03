#########################################
#
# Shiny Web Application: Interactive Network Visualization
#
# Purpose: 
# This application visualizes pathway enrichment similarity networks
# from multi-strain/multi-condition analysis. Users can interactively
# filter networks based on conservation across strains and customize
# visual appearance.
#
# Input files:
# - Node file (.tsv): contains pathway information, enrichment statistics
# - Edge file (.tsv): contains pathway-pathway similarity scores
#
# Key features:
# - Filter networks by strain conservation
# - Customize color schemes
# - Interactive visualization with ggraph
#
# How to run:
# 1. Click 'Run App' button in RStudio
# 2. Upload node and edge .tsv files
# 3. Adjust filtering and coloring parameters
# 4. Click 'Plot Network' to visualize
#
# Reference: https://shiny.rstudio.com/
#
#########################################

# Load required packages
library(shiny)          # Web application framework
library(shinythemes)    # Bootstrap themes for Shiny
library(igraph)         # Network analysis and manipulation
library(ggraph)         # Grammar of graphics for networks
library(tidyverse)      # Data manipulation (dplyr, ggplot2, etc.)
library(tidygraph)      # Tidy interface to igraph

#########################################
# USER INTERFACE (UI) DEFINITION
#########################################
# The UI defines what users see and interact with

ui <- fluidPage(
  theme = shinytheme("cerulean"),  # Apply cerulean Bootstrap theme
  titlePanel("Network Graph from Node and Edge Files"),
  
  # Create sidebar layout (sidebar + main panel)
  sidebarLayout(
    
    ## SIDEBAR PANEL: Input controls
    sidebarPanel(
      
      # File upload: Node table
      # Nodes represent pathways with attributes (name, p-value, strain count, etc.)
      fileInput("nodeFile", 
                "Upload Node Table (.tsv)", 
                accept = ".tsv"),
      
      # File upload: Edge table
      # Edges represent similarity between pathways
      fileInput("edgeFile", 
                "Upload Edge Table (.tsv)", 
                accept = ".tsv"),
      
      # Numeric input: Filter by number of strains
      # Shows pathways conserved across N strains/conditions
      numericInput("numStrains", 
                   "Number of strains/conditions to show:", 
                   value = 5, 
                   min = 1, 
                   max = 9, 
                   step = 1),
      
      # Radio buttons: Choose filter type
      # "Exact": show pathways in exactly N strains
      # "Less Than or Equal": show pathways in ≤ N strains
      radioButtons("filterType", 
                   "Filter Type (#strains/#conditions):", 
                   choices = c("Exact" = "exact", 
                               "Less Than or Equal" = "lte"),
                   selected = "lte"),
      
      # Radio buttons: Choose color scheme type
      radioButtons("colorChoice", 
                   "Color Palette Option:",
                   choices = c("Preset Palette" = "preset", 
                               "Custom Hex Codes" = "custom")),
      
      # Conditional panel: Show preset palette options
      # Only visible when "Preset Palette" is selected
      conditionalPanel(
        condition = "input.colorChoice == 'preset'",
        selectInput("paletteName", 
                    "Choose Preset Palette:",
                    choices = c("Blues", "Reds", "Viridis", "Plasma"))
      ),
      
      # Conditional panel: Custom hex color codes
      # Only visible when "Custom Hex Codes" is selected
      conditionalPanel(
        condition = "input.colorChoice == 'custom'",
        textInput("hexColors", 
                  "Enter Hex Codes (comma-separated):",
                  "#4575B4,#91BFDB,#FEE090,#FC8D59,#D73027")
      ),
      
      # Action button: Trigger network plotting
      # Network only updates when this button is clicked
      actionButton("plotButton", "Plot Network")
    ),
    
    ## MAIN PANEL: Display output
    mainPanel(
      # Plot output: Displays the network graph
      plotOutput("networkPlot", 
                 width = "100%", 
                 height = "800px")
    )
  )
)

#########################################
# SERVER LOGIC DEFINITION
#########################################
# The server defines how the app responds to user inputs

server <- function(input, output) {
  
  ## REACTIVE: Load node data from uploaded file
  # Executes whenever user uploads a new node file
  # Returns data frame with pathway information
  nodeData <- reactive({
    req(input$nodeFile)  # Require file to be uploaded
    read.delim(input$nodeFile$datapath, 
               header = TRUE, 
               sep = "\t", 
               stringsAsFactors = FALSE)
  })
  
  ## REACTIVE: Load edge data from uploaded file
  # Executes whenever user uploads a new edge file
  # Returns data frame with pathway similarity scores
  edgeData <- reactive({
    req(input$edgeFile)  # Require file to be uploaded
    read.delim(input$edgeFile$datapath, 
               header = TRUE, 
               sep = "\t", 
               stringsAsFactors = FALSE)
  })
  
  ## EVENT REACTIVE: Process data when "Plot Network" is clicked
  # Only updates when user clicks the plot button
  # Filters nodes and edges based on user selections
  graphData <- eventReactive(input$plotButton, {
    
    # Get uploaded data
    nodes <- nodeData()
    edges <- edgeData()
    num <- input$numStrains
  
    # Filter nodes based on strain conservation
    # This controls which pathways appear in the network
    if(input$filterType == "exact") {
      # Show pathways present in exactly N strains
      nodes_filtered <- nodes %>% 
        filter(nb_strain_with_node == num)
    } else {
      # Show pathways present in ≤ N strains (more inclusive)
      nodes_filtered <- nodes %>% 
        filter(nb_strain_with_node <= num)
    }
   
    # Filter edges: only keep connections between retained nodes
    # This ensures network integrity (no orphan edges)
    edge_nodes <- nodes_filtered$name
    edges_filtered <- edges %>%
      filter(from %in% edge_nodes & to %in% edge_nodes)
    
    # Create igraph network object
    # Nodes = pathways, Edges = similarity relationships
    # directed=FALSE: similarity is symmetric (A~B equals B~A)
    ig <- igraph::graph_from_data_frame(
      d = edges_filtered,          # Edge data frame  
      vertices = nodes_filtered,   # Node data frame
      directed = FALSE
    )
    
    # Return list containing graph and data
    return(list(graph = ig, 
                nodes = nodes_filtered, 
                edges = edges_filtered))
  })
  
  ## OUTPUT: Render network plot
  # Creates the visual network representation
  output$networkPlot <- renderPlot({
    
    # Require graph data to be ready
    req(graphData())
    
    # Extract components from graphData
    ig <- graphData()$graph
    nodes <- graphData()$nodes
    edges <- graphData()$edges
    
    # Convert igraph to tidygraph for better integration with ggplot2
    # Add node labels for visualization
    set.seed(1234)  # Reproducible layout
    tg <- tidygraph::as_tbl_graph(ig) %>% 
      tidygraph::activate(nodes) %>% 
      dplyr::mutate(label = name)  # Add label column for node names
    
    ## Prepare color palette based on user selection
    if (input$colorChoice == "preset") {
      # Use preset color palettes from RColorBrewer or viridis
      palette_choice <- switch(
        input$paletteName,
        "Blues" = RColorBrewer::brewer.pal(9, "Blues"),
        "Reds" = RColorBrewer::brewer.pal(9, "Reds"),
        "Viridis" = viridis::viridis(9),
        "Plasma" = viridis::plasma(9)
      )
      # Interpolate to match number of strains
      color_palette <- colorRampPalette(palette_choice)(input$numStrains)
    } else {
      # Parse custom hex codes from user input
      # Split comma-separated string and trim whitespace
      color_palette <- strsplit(input$hexColors, ",")[[1]]
      color_palette <- trimws(color_palette)  
    }
    
    ## Create network visualization using ggraph
    # ggraph extends ggplot2 for network data
    # Reference: https://ggraph.data-imaginist.com/
    plot_tg <- tg %>%
      ggraph::ggraph(layout = "fr") +  # Fruchterman-Reingold force-directed layout
      
      # Draw edges (connections between pathways)
      # Arc shape emphasizes network structure
      geom_edge_arc(
        colour = "gray50",              # Gray edges
        lineend = "round",              # Rounded line ends
        strength = 0.1,                 # Arc curvature
        aes(edge_width = edges$nb_strain_with_edge,    # Edge thickness by conservation
            alpha = edges$nb_strain_with_edge)         # Edge transparency by conservation
      ) +
      
      # Draw node points (pathways)
      # Size represents strain conservation
      # Color represents strain count (gradient)
      ggraph::geom_node_point(
        size = nodes$nb_strain_with_node * 1.5,  # Larger nodes = more conserved
        aes(color = nodes$nb_strain_with_node)   # Color gradient by conservation
      ) +
      
      # Add node border (second layer with hollow shape)
      # Creates visual depth and improves visibility
      ggraph::geom_node_point(
        size = nodes$nb_strain_with_node * 1.5, 
        shape = 21  # Hollow circle
      ) + 
      
      # Add text labels for pathway names
      # ggrepel prevents label overlap
      geom_node_text(
        aes(label = name), 
        repel = TRUE,                          # Avoid overlapping labels
        point.padding = unit(0.2, "lines"),    # Padding around points
        colour = "gray10"                      # Dark gray text
      ) +
      
      # Apply color gradient to nodes
      # Maps strain count to color scale
      scale_color_gradientn(
        colors = color_palette,
        limits = c(1, input$numStrains),
        breaks = 1:input$numStrains
      ) + 
      
      # Customize edge width range
      scale_edge_width(range = c(0.5, 2.5)) +
      
      # Customize edge transparency range
      scale_edge_alpha(range = c(0.2, 0.5)) +
      
      # Use clean white background
      theme_graph(background = "white") +
      
      # Position legend at bottom
      theme(
        legend.position = "bottom", 
        legend.title = element_blank()  # Remove legend title
      ) 
      
    # Return the plot
    plot_tg
    
  })  # End renderPlot
  
}  # End server function

#########################################
# RUN THE APPLICATION
#########################################
# This line launches the Shiny app
# Combines UI and server into a functional web application
shinyApp(ui = ui, server = server)
