
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

rm(list=ls())
Sys.getenv('HOME')
Sys.getenv('BIOCONDUCTOR_CONFIG_FILE')
Sys.setenv(BIOCONDUCTOR_CONFIG_FILE="http://bioconductor.org/config.yaml")
chooseBiocMirror()
BiocManager::install(c("graph", "RBGL", "combinat", force=TRUE))
install.packages('DESeq2')
library(DESeq2)
# find packages
BiocManager::available()

#install.packages('https://cran.r-project.org/src/contrib/Archive/QuACN/QuACN_1.8.0.tar.gz', contriburl=NULL, type="source")
install.packages('rlang')
install.packages("Rpmf", dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages("https://cran.r-project.org/src/contrib/Archive/QuACN/QuACN_1.8.0.tar.gz", repos = NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Rmpfr_0.9-2.tar.gz", repos = NULL, type="source")
install.packages('https://www.bioconductor.org/packages/release/bioc/src/contrib/graph_1.78.0.tar.gz', repos=NULL, type="source")
install.packages('https://cran.r-project.org/src/contrib/combinat_0.0-8.tar.gz', repos=NULL, type="source")
install.packages('https://cran.r-project.org/src/contrib/gmp_0.7-1.tar.gz')
#install.packages('https://cran.r-project.org/web/packages/multiplex/index.html#:~:text=multiplex_3.1.0.tar.gz')
install.packages('multiplex')

install.packages('Rgraphviz')
library(Rgraphviz)
remove.packages(c("BiocVersion", "BiocInstaller")) # repeat 'till all removed
install.packages("BiocManager")
if(!require('QuACN')) {
  BiocManager::install('QuACN')
  library('QuACN')
}



setRepositories() #Bioconductor
BiocManager::install.packages("graph")
install.packages("RBGL")
install.packages("combinat")
install.packages('QuACN')
install.packages("graphNEL")
library(RBGL)
library(combinat)
library(graph)
library(igraph)
library(QuACN)
# for edgeL
library(multiplex)

set.seed(123)
set.seed(666)
g <- randomGraph(1:8, 1:5, 0.36, weights=FALSE)


# degreeCentrality example
#edges <- c(1,2,1,4,2,3,2,4,3,4,4,5,5,6,6,7,6,8,6,9)
#g<-graph(edges, 9,directed=FALSE) 
#l<-matrix(c(1,1,0,0,1,-1,2,0,3,0,4,0,5,-1,5,0,5,1),nrow = 9,ncol=2,byrow=TRUE) #custom layout
#plot(g,layout = l)
#par(mfrow=c(1,2))
#plot(g,layout=l)
#plot(g,vertex.label=round(degree(g)/sum(degree(g)),2),layout=l,vertex.label.dist=2.5,vertex.size=80*degree_centrality)

#degree_centrality<-degree(g)/sum(degree(g))
#degree_centrality

#closeness(g)
#par(mfrow=c(1,2))
#plot(g,layout=l)
#plot(g,vertex.label=round(closeness(g),2),layout=l,vertex.label.dist=2.5,vertex.size=190*closeness(g))


# read .edges Graph from networkrepository.com
# NOTE graphNEL are typically used for directed graphs
listOfEdgesFiles <- list("bio-CE-GT", "bio-CE-HT", "bio-CE-LC", "bio-DM-HT", "bio-grid-mouse", "bio-grid-plant", "bio-yeast", "bio-yeast-protein-inter")

# read graph and simplify the iGraph
wikiGraph <- read.graph("/Users/eldoradoyankee/Library/CloudStorage/OneDrive-FFHS/Bachelorthesis/BT/Graphen/ENZYMES_g224/ENZYMES_g224.edges", format=c("ncol"))
# https://stackoverflow.com/questions/25899516/importing-edge-list-in-igraph-in-r
graph <- as.directed(read.graph("/Users/eldoradoyankee/Library/CloudStorage/OneDrive-FFHS/Bachelorthesis/BT/Graphen/biological/bio-CE-HT/bio-CE-HT.edges", format = c("ncol"), directed=TRUE, weights="yes"), mode=c("directed"))
graph <- read.graph("/Users/eldoradoyankee/Library/CloudStorage/OneDrive-FFHS/Bachelorthesis/BT/Graphen/biological/bio-grid-plant/bio-grid-plant.edges", format = c("ncol"))
dd <- read.table("/Users/eldoradoyankee/Library/CloudStorage/OneDrive-FFHS/Bachelorthesis/BT/Graphen/biological/bio-CE-GT/bio-CE-GT.edges", header=FALSE)
#openAI - 
graph_obj <- new("graphNEL", nodes = character(0), edgemode = "undirected")
graph_obj <- graph(edges=NULL)
graph_obj <- graphNEL()
# where all edges are stored for is_edge_added function
added_edges <- list()

# Loop through the edges data and add the nodes and edges to the graphNEL object

dd <- matrix(c(1, 2, 1, 1, 3, 1, 1, 4, 1), ncol=3, byrow=TRUE)
print(dd)

for (i in 1:nrow(dd)) {
  node1 <- as.character(dd[i, 1])
  node2 <- as.character(dd[i, 2])
  weight <- as.numeric(dd[i, 3])
  
  
  # Add nodes to the graph if they don't already exist
  #if (!(node1 %in% V(graph_obj))) {
  #  graph_obj <- add_vertices(graph_obj, node1)
  #}
  #if (!(node2 %in% V(graph_obj))) {
  #  graph_obj <- add_vertices(graph_obj, node2)
  #}
  if (!(node1 %in% nodes(graph_obj))) {
    graph_obj <- addNode(node1, graph_obj)
  }
  if (!(node2 %in% nodes(graph_obj))) {
    graph_obj <- addNode(node2, graph_obj)
  }
  
  
  # Add the edge to the graph
  #attr = list(name = c(node1, node2))
  
  #if (!edgeExists(node1, node2, graph_obj) && !edgeExists(node2, node1, graph_obj)) {
  #  graphNEL_obj <- addEdge(node1, node2, edgeData = weight, graph_obj)
  #  graphNEL_obj <- addEdge(node2, node1, edgeData = weight, graph_obj)
  #}
  
  if (!is_edge_added(node1, node2, added_edges) && !is_edge_added(node2, node1, added_edges)) {
    graph_obj <- addEdge(graph_obj, from = node1, to = node2)
    graph_obj <- addEdge(graph_obj, from = node2, to = node1)
    added_edges <- c(added_edges, list(edge(node1, node2)))
  }
  
  
  #graph_obj <- add_edges(graph_obj, c(node1, node2), weight=weight)
  #graph_obj <- add_edges(graph_obj, c(node2, node1), weight=weight)
  
  #graph_obj <- graph_obj + edge(node1, node2, weight = weight)
  #E(graph_obj)$weight <- weight
  #graph_obj <- add_edges(graph_obj, node1, node2, weight=weight, .names = "weight")
}
print(graph_obj)
graph_obj <- getLargestSubgraph(graph_obj)
isFullyConnected(graph_obj)

is_fully_connected <- function(graph_obj) {
  return(is.connected(graph_obj))
}
fully_connected <- is_fully_connected(graph_obj)
#graph <- graph_from_graphNEL(graph_obj)
print(is.igraph(graph_obj))
graph_obj <- new("graphNEL", nodes=as.character(vertices(graph_obj)), edgeL=list(edges(graph_obj)), edgemode = "undirected")


print(graph_obj)
plot(graph_obj)
#graph_obj <- as_graphNEL(graph_obj)
wiener(graph_obj)
###
###

wikiGraph <- simplify(wikiGraph, remove.multiple = TRUE, remove.loops = TRUE)
plot(wikiGraph)

validateGraph(graph)


g <- as(A, "graphNEL")
g <- randomGraph(1:8, 1:5, 0.36, weights=FALSE)

wiener(g)


#https://stackoverflow.com/questions/25721884/how-should-i-deal-with-package-xxx-is-not-available-for-r-version-x-y-z-wa
sessionInfo()


# check if an edge is already added
is_edge_added <- function(node1, node2, added_edges) {
  for (edge in added_edges) {
    if (identical(edge, edge(node1, node2))) {
      return(TRUE)
    }
  }
  return(FALSE)
}


dfs <- function(graph, start) {
  visited <- rep(FALSE, length(nodes(graph)))
  visited[start] <- TRUE
  
  stack <- start
  while (length(stack) > 0) {
    current <- stack[length(stack)]
    stack <- stack[-length(stack)]
    neighbors <- neighbors(graph, nodes(graph)[current])
    for (neighbor in neighbors) {
      if (!visited[neighbor]) {
        visited[neighbor] <- TRUE
        stack <- c(stack, neighbor)
      }
    }
  }
  all(visited)
}

# Test if graph is fully connected
isFullyConnected <- function(graph) {
  nodes <- nodes(graph)
  if (length(nodes) <= 1) {
    return(TRUE)
  }
  start_node <- nodes[1]
  dfs_result <- dfs(graph, root=start_node)
  all(nodes %in% dfs_result)
}
