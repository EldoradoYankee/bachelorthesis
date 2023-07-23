
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
install.packages("gRbase")

install.packages('Rgraphviz')
library(Rgraphviz)
remove.packages(c("BiocVersion", "BiocInstaller")) # repeat 'till all removed
install.packages("BiocManager")
if(!require('QuACN')) {
  BiocManager::install('QuACN')
  library('QuACN')
}

#https://stackoverflow.com/questions/25721884/how-should-i-deal-with-package-xxx-is-not-available-for-r-version-x-y-z-wa
sessionInfo()

setRepositories() #Bioconductor
BiocManager::install.packages("graph")
install.packages("RBGL")
install.packages("combinat")
install.packages('QuACN')
install.packages("graphNEL")
install.packages("Seurat")
library(RBGL)
library(combinat)
library(graph)
library(igraph)
library(QuACN)
library(gRbase)
library(Seurat)
library(Matrix)
# for edgeL
library(multiplex)

# read .edges Graph from networkrepository.com
# NOTE graphNEL are typically used for directed graphs
# "bio-CE-GT", "bio-CE-HT", "bio-CE-LC", "bio-DM-HT","bio-grid-mouse", "bio-grid-plant", "bio-yeast-protein-inter"
# "fb-pages-company", "fb-pages-government", "fb-pages-politician", "fb-pages-sport", fb-pages-tvshows", "soc-advogato", "soc-anybeat", "soc-hamsterster"
listOfEdgesFiles <- list( "bio-CE-GT", "bio-CE-HT", "bio-CE-LC", "bio-DM-HT", "bio-grid-mouse", "bio-grid-plant", "bio-yeast-protein-inter")
listOfEdgesFiles <- list("eins")

# list to fill with all the graphs and their results from the applied graph measures 
wienerIndices <- list()
df <- data.frame(matrix(ncol = length(listOfEdgesFiles), nrow = 0))
colnames(df) <- listOfEdgesFiles
print(df)

###
# start loop for wienerIndex with multiple graphs
for (edgeFile in listOfEdgesFiles) {
  # set path
  pathVar = paste("/Users/eldoradoyankee/Library/CloudStorage/OneDrive-FFHS/Bachelorthesis/BT/Graphen/biological/", edgeFile, "/", edgeFile, ".edges", sep="")
  print(pathVar)
  
  
  # check graph for csv
  #is_csv <- tryCatch({
  #  read.ftable(grepl(",", pathVar))
  #  print("table read as csv")
  #  TRUE
  #}, error = function(e) {
  #  FALSE
  #})
  #print(paste("the dataset is csv:", is_csv))
  #read graph
  dd <- read.table(pathVar, header = FALSE)
  #print(dd)
  print("table read")
  
  
  # read graph
  #if (!is_csv) {
  #  dd <- read.table(pathVar, header=FALSE)
  #  print("datatable non csv")
  #} else {
  #  dd <- read.table(pathVar, header=FALSE, sep = ",", row.names = NULL)
  #  print("datatable csv")
  #}
  
  # example graph for testing
  #dd <- matrix(c(1, 2, 1, 1, 3, 1, 1, 4, 5), ncol=3, byrow=TRUE)

  # reset every time for new graph measure or graph
  graph_obj <- graphNEL(edgeL=list(), edgemode = "undirected")
  
  # where all edges are stored for is_edge_added function and check if edge already exists between vertices - also reset everytime
  added_edges <- list()

  #print(paste("length of the edgeList", length(listOfEdgesFiles)))
  #edL2 <- vector("list", length=length(listOfEdgesFiles))
  #print(paste("this is the edgelist", edL2))
  
  # Loop through the edges data and add the nodes and edges to the graphNEL object
  for (i in 1:nrow(dd)) {
    node1 <- as.character(dd[i, 1])
    node2 <- as.character(dd[i, 2])
    #print(dd[i, 3]) #prints the weight found in 3rd column @ position i
    if (ncol(df) == 3 && is.null(dd[i, 3])) {
      weight <- as.numeric(dd[i, 3])
    }
    
    if (!(node1 %in% nodes(graph_obj))) {
      graph_obj <- addNode(node1, graph_obj)
    }
    if (!(node2 %in% nodes(graph_obj))) {
      graph_obj <- addNode(node2, graph_obj)
    }
    
    
    # Add the edge to the graph
    if (!is_edge_added(node1, node2, added_edges) && !is_edge_added(node2, node1, added_edges)) {
      graph_obj <- addEdge(graph_obj, from = node1, to = node2)
      graph_obj <- addEdge(graph_obj, from = node2, to = node1)
      added_edges <- c(added_edges, list(edge(node1, node2)))
      #print("edge added to list")
      #edgeData(graph_obj, node1, node2)$weight <- weight
      #print("weight added")
      
    }
    
    
    # add edgeList to graphNEL object
    #edL2[[i]] <- list(edges=c(node1, node2)[i], weights=weight)
    #print(paste("length of edL2", length(edL2)))
    #print(paste("the edL2 edgeList", edL2))
  }
  #print(paste("all the added weights", edL2))
  
  
  
  
  # Set the weight attribute of the edge
  #edge_obj <- edge(graph_obj, from = node1, to = node2)
  #edgeDataDefaults(graph_obj, attr = "weight")
  #setEdgeData(graph_obj, edge_obj, attr = "weight", value = weight)
  
  #graph_obj@edgeL <- edL2
  #print("graph_obj@edgeL edl2 added")
  # Set the weight attribute of the edge
  #edge_obj <- edge(graph_obj, from = node1, to = node2)
  #edgeDataDefaults(graph_obj, attr = "weight")
  #setEdgeData(graph_obj, edge_obj, attr = "weight", value = weight)

  
  #if (is.connected(graph_from_graphnel(graph_obj))) {
    graph_obj <- getLargestSubgraph(graph_obj)
  #}
  # not all graphs are connected
  #if(!isFullyConnected(graph_obj)) {
  #  graph_obj <- getLargestSubgraph(graph_obj)
  #}
  
  wienerIndices[[length(wienerIndices)+1]] <- wiener(graph_obj)
  print(paste("The WienerIndex of the current Graph", edgeFile, "is", wiener(graph_obj)))
  
  plot(graph_obj)
  #print(edgeL(graph_obj))
}
df[nrow(df) + 1,] = c(wienerIndices)


print(wienerIndices)
print(df)

isFullyConnected(graph_obj)

# is fully Connected graphNEL??
is_fully_connected <- function(graph_obj) {
  return(is.connected(graph_obj))
}
fully_connected <- is_fully_connected(graph_obj)
print(fully_connected)

# get and print weights from graphNEL object
edges <- edges(graph_obj)
for (edge in edges) {
  print(edgeWeights(edge))
}

print(graph_obj)
plot(graph_obj)
###

wikiGraph <- simplify(wikiGraph, remove.multiple = TRUE, remove.loops = TRUE)
plot(wikiGraph)

validateGraph(graph)








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
isFullyConnected <- function(graph_obj) {
    num_nodes <- vcount(graph_obj)
    
    if (num_nodes <= 1) {
      return(TRUE)
    }
    
    for (i in 1:num_nodes) {
      for (j in (i + 1):num_nodes) {
        if (!is.connected(graph_obj, nodes = c(i, j))) {
          return(FALSE)
        }
      }
    }
    
    return(TRUE)
  }