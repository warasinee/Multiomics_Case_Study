###############################
#
# Network visualization using Cytoscape - DIABLO results 
#
# Written by Warasinee Mujchariyakul
#
###############################

# Load R packages
library(RCy3)
library(igraph)
set.seed(123) # for reproducibility

###############################

## Introduction to Cytoscape in R ##
# Ref: http://127.0.0.1:13113/session/Rvig.1006a4f648dea.html
# Overview of RCy3 : http://127.0.0.1:13113/library/RCy3/doc/Overview-of-RCy3.html

library(RCy3)

# Open cytoscape software!!!!!

# Getting started (Connecting R to Cytoscape)
cytoscapePing () 

# Check version (Here: apiVersion (V1), cytoscapeVersion (3.9.1))
cytoscapeVersionInfo ()

# More information about package RCy3
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