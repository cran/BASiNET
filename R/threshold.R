#'@title Applies threshold on the network from a value
#'@name threshold
#'
#'@description Given an integer value X, a cut, that is, edges that are cut will be assigned zero. This cut will be done in the network where the edges have a weight less than the value of X.
#'
#'@param x Value that would limit the edges
#'@param net Network where the edges will be cut
#'
#'
#'@return Returns the complex network with the cuts already made
#'
#' @author Eric Augusto Ito
#'
#'
#'
#' @import igraph

threshold <- function(x, net){
	# tamanho <- length(vetor<-net[1,])
	# for(i in 1:tamanho){
	# 	for(k in 1:tamanho){
	# 		if((net[i,k]!=0)&&(net[i,k]<x)){
	# 			net[i,k]<-0		
	# 		}
	# 	}
	# }

	if(x==1){

	}else{
		matriz<-as_adjacency_matrix(net)
		matriz[(matriz==(x-1))]<-0
		net<-graph_from_adjacency_matrix(matriz, mode="undirected")
	}
	return(net)
}