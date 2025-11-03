###############################
#
# Network visualization using Cytoscape - DIABLO results 
#
# This script demonstrates how to:
# 1. Connect R to Cytoscape using RCy3
# 2. Import DIABLO network (GML file from mixOmics)
# 3. Customize network appearance (nodes, edges, layout)
# 4. Export publication-ready network figures
#
# Prerequisites:
# - Cytoscape software must be installed and running
# - GML file generated from DIABLO analysis
#
# Reference: http://mixomics.org/mixomics-study-examples/diablo-tcga-case-study/#visualize-in-cytoscape
# RCy3 Documentation: https://cytoscape.org/RCy3/
#
# Written by Warasinee Mujchariyakul
#
###############################

# Load R packages
# RCy3: R interface to Cytoscape for programmatic network visualization
library(RCy3)

# igraph: network analysis and file I/O for GML format
library(igraph)

# Set seed for reproducibility in layout algorithms
set.seed(123)

###############################

## STEP 1: Connect R to Cytoscape ##
# Before running this script, make sure Cytoscape is open and running!
# Download Cytoscape: https://cytoscape.org/download.html
# Reference: https://cytoscape.org/RCy3/articles/Overview-of-RCy3.html

library(RCy3)

# Test connection to Cytoscape
# This should return "You are connected to Cytoscape!" if successful
# If it fails, ensure Cytoscape is running and port 1234 is not blocked
cytoscapePing() 

# Check Cytoscape version information
# Verifies API version and Cytoscape version for compatibility
# Example output: apiVersion (v1), cytoscapeVersion (3.9.1)
cytoscapeVersionInfo()

# Access RCy3 package documentation and tutorials
# Opens vignettes in your browser with detailed examples
browseVignettes("RCy3")

###############################

# Import .gml file (output from MixOmics (DIABLO))

library(igraph)
my.igraph.network <- igraph::read_graph("/Users/wmujchariyak/Desktop/myNetwork_conserved_Strep_05.gml",format=c("gml"))

###############################

# From igraph to Cytoscape (Now you will see the network in Cytoscape)
createNetworkFromIgraph(my.igraph.network,"myIgraph.Strep")

# Verify the number of nodes (26) and edges (80):
igraph::vcount(my.igraph.network)
igraph::ecount(my.igraph.network)

# Common iGraph functions
# Calculate degree for all nodes
degAll <- igraph::degree(my.igraph.network, v = igraph::V(my.igraph.network), mode = "all")

# Let’s decide on a layout
getLayoutNames()
#layoutNetwork('force-directed defaultSpringLength=400 defaultSpringCoefficient=0.000003')
layoutNetwork('hierarchical')
###############################
## Adjust Node (Nodes = Selected variables from DIABLO) ##

# Change node shape
getNodeShapes()
setNodeShapeDefault("ELLIPSE")

# Get the names of all nodes 
getAllNodes()

# Get the names of all columns in a table
getTableColumnNames(table = "node")

# Set node's color based on data in "group" column
# For more color Hex code : https://www.rapidtables.com/web/color/RGB_Color.html
setNodeColorMapping(table.column = "group" , table.column.values = c("Transcript", "Protein", "Metabolite_GC"), mapping.type = "discrete", colors=c('#E5CCFF','#FFCCCC','#CCFFCC'))

###############################
## Adjust Edge (Edges = Interaction between variables) ##

# Get the names of all edges
getAllEdges()

# Get the names of all columns in a table
getTableColumnNames(table = "edge")
getEdgeCount()

# Set edge's color based on correlation data in "weight" column
setEdgeColorMapping(table.column = "weight", c(-1,1), c('#000000','#FF9933'))

###############################
# Adjust network 

# Calculate degree for all nodes
degree_values <- igraph::degree(my.igraph.network, v = igraph::V(my.igraph.network), mode = "all")

# Loads data (e.g., degree_values) into Cytoscape node tables
df_degree <- data.frame(degree_values)
loadTableData(df_degree, data.key.column = "row.names", table = "node", table.key.column = "name")

# Set size & font of node
#setNodeSizeMapping(table.column = "degree_values",sizes = c(150,270), mapping.type = "c")
#setNodeFontSizeMapping(table.column = "degree_values",sizes = c(60,65), mapping.type = "c")
setNodeSizeMapping(table.column = "degree_values",sizes = c(80,150), mapping.type = "c")
setNodeFontSizeMapping(table.column = "degree_values",sizes = c(30,45), mapping.type = "c")

# Note: You can manually adjust network in Cytoscape software directly
###############################
# Save 

# Cytoscape format (.cys) 
saveSession("/Users/wmujchariyak/Desktop/my.network.Strep.081025.cys") 

# PNG/PDF/SVG format
?exportImage
full.path=paste("/Users/wmujchariyak/Desktop/","my.network.Strep.081025",sep="/")
exportImage(full.path, "PNG", zoom=200) #.png scaled by 200%
exportImage(full.path, "PDF") #.pdf
exportImage(full.path, "SVG") 

###############################
# Note: Manually adjust network (SVG format) using Inkscape software  