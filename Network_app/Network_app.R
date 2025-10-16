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

# UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                titlePanel("Network Graph from Node and Edge Files"),
                sidebarLayout(
                  sidebarPanel(
                    fileInput("nodeFile", "Upload Node Table (.tsv)", accept = ".tsv"),
                    fileInput("edgeFile", "Upload Edge Table (.tsv)", accept = ".tsv"),
                    numericInput("numStrains", "Number of strains/conditions to show:", 5, min=1, max=9, step=1),
                    
                    # Filter type: exact or less than/equal
                    radioButtons("filterType", "Filter Type (#strains/#conditions):", 
                                 choices = c("Exact" = "exact", "Less Than or Equal" = "lte"),
                                 selected = "lte"),
                    
                    # Choose color palette or hex codes
                    radioButtons("colorChoice", "Color Palette Option:",
                                 choices = c("Preset Palette" = "preset", "Custom Hex Codes" = "custom")),
                    # If preset, show selectInput
                    conditionalPanel(
                      condition = "input.colorChoice == 'preset'",
                      selectInput("paletteName", "Choose Preset Palette:",
                                  choices = c("Blues", "Reds", "Viridis", "Plasma"))
                    ),
                    # If custom, show textInput
                    conditionalPanel(
                      condition = "input.colorChoice == 'custom'",
                      textInput("hexColors", "Enter Hex Codes (comma-separated):",
                                "#4575B4,#91BFDB,#FEE090,#FC8D59,#D73027")
                    ),
                    
                    actionButton("plotButton", "Plot Network")
                  ),
                  
                  mainPanel(
                    plotOutput("networkPlot")
                  )
                )
)

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
    num <- input$numStrains
  
    # Filter nodes based on filter type
    if(input$filterType == "exact") {
      nodes_filtered <- nodes %>% filter(nb_strain_with_node == num)
    } else {
      nodes_filtered <- nodes %>% filter(nb_strain_with_node <= num)
    }
   
    # Filter edges that connect only to these nodes
    edge_nodes <- nodes_filtered$name
    edges_filtered <- edges %>%
      filter(from %in% edge_nodes & to %in% edge_nodes)
    
    ig <- igraph::graph_from_data_frame(d=edges_filtered, vertices=nodes_filtered, directed=FALSE)
    return(list(graph=ig, nodes=nodes_filtered, edges=edges_filtered))
    
  })
  # Plot using ggraph
  output$networkPlot <- renderPlot({
    req(graphData())
    ig <- graphData()$graph
    nodes <- graphData()$nodes
    edges <- graphData()$edges
    
    
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
      color_palette <- colorRampPalette(palette_choice)(input$numStrains)
    } else {
      # Custom hex codes → split by comma
      color_palette <- strsplit(input$hexColors, ",")[[1]]
      color_palette <- trimws(color_palette)  
    }
    
    
  
    # 4. Plot the network
    plot_tg <- tg %>%
      ggraph::ggraph(layout = "fr") +
      # Add arcs for edges with various aesthetics
      geom_edge_arc(colour = "gray50",
                    lineend = "round",
                    strength = .1,
                    aes(edge_width = edges$nb_strain_with_edge,
                        alpha = edges$nb_strain_with_edge)) +
      # Add points for nodes with size based on nb_strain_with_node and color based on gene_set_size
      ggraph::geom_node_point(size = nodes$nb_strain_with_node*1.5,
                              aes(color = nodes$nb_strain_with_node)) +
      
      # Add edge of nodes on the next layer (shape = 21)
      ggraph::geom_node_point(size = nodes$nb_strain_with_node*1.5, shape = 21) + 
      
      # Add text labels for nodes with various aesthetics
      geom_node_text(aes(label = name), 
                     repel = TRUE, 
                     point.padding = unit(0.2, "lines"), 
                     colour = "gray10") +
      
      # Color  
      scale_color_gradientn(colors = color_palette,
                            limits = c(1, input$numStrains),
                            breaks = 1:input$numStrains) + 
      
      # Adjust edge width and alpha scales
      scale_edge_width(range = c(0.5, 2.5)) +
      scale_edge_alpha(range = c(0.2, 0.5)) +
      
      # Set graph background color to white
      theme_graph(background = "white") +
      
      # Adjust legend position to the top
      theme(legend.position = "bottom", 
            # suppress legend title
            legend.title = element_blank()) 
      
    plot_tg
    
  })
}
# Run the app
shinyApp(ui = ui, server = server)
