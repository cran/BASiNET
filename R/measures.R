#'@title Abstracting Characteristics on Network Structure
#'@name measures
#'
#'@description Given a graph, it is made up of several features on the graph structure and returns a vector with the data obtained
#'
#'@param graph The complex network that will be measured
#'
#'
#'@return Return a vector with the results of the measurements in order:
#'Average shortest path length, clustering Coefficient, degree, assortativity,
#'betweenness, standard deviation, maximum, minimum, number of motifs
#'size 3 and number of motifs of size 4
#'
#' @author Eric Augusto Ito
#'
#' 
#' @import igraph
#' @importFrom stats sd

measures <- function(graph){
	measures<-c()
	measures<-c(measures,average.path.length(graph,directed=FALSE,unconnected=FALSE))
	measures<-c(measures,transitivity(graph,type=c("undirected"),vids=NULL,weights=NULL, isolates=c("NaN","zero")))
	measures<-c(measures,mean(degree(graph, v=V(graph), normalized=FALSE)))
	measures<-c(measures,assortativity_degree(graph, directed = FALSE))
	measures<-c(measures,mean(betweenness(graph, v = V(graph), directed = FALSE, weights = NULL, normalized = FALSE)))
	measures<-c(measures,sd(degree(graph, v=V(graph), normalized= FALSE), na.rm = FALSE))
	measures<-c(measures,which.max(degree(graph, v=V(graph), normalized=FALSE)))
	measures<-c(measures,which.min(degree(graph, v=V(graph), normalized=FALSE)))
	measures<-c(measures,(count_motifs(graph, size = 3)))
	measures<-c(measures,(count_motifs(graph, size = 4)))
	names(measures)[1]<-"ASPL"
	names(measures)[2]<-"CC"
	names(measures)[3]<-"DEG"
	names(measures)[4]<-"ASS"
	names(measures)[5]<-"BET"
	names(measures)[6]<-"SD"
	names(measures)[7]<-"MAX"
	names(measures)[8]<-"MIN"
	names(measures)[9]<-"MT3"
	names(measures)[10]<-"MT4"

	return(measures)
}
