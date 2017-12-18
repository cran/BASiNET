#'@title Creates an undirected graph from a biological sequence
#'@name createGraph2D
#'
#'@description A function that from a biological sequence generates a graph not directed having as vertices words, this being able to have its size parameter set by the 'word' parameter. The connections between words depend of the 'step' parameter that indicates the next connection to be formed.
#'
#'@param matrix matrix of the measure for the creation of two-dimensional graph
#'@param numSeqMRNA number of mRNA sequences
#'@param numSeqLNCRNA number of lncRNA sequences
#'@param nameMeasure name of the measure to put in the title of the graph
#'
#'
#'@return Returns the non-directed graph formed through the sequence
#'
#' @author Eric Augusto Ito
#'
#'
#' @import igraph
#' @importFrom grDevices dev.new
#' @importFrom graphics lines plot title


createGraph2D <- function(matrix, numSeqMRNA,numSeqLNCRNA, nameMeasure){

	dev.new()
	somador<-(1/(length(matrix[1,])-1))
	threshold<-seq(0,1,somador)
	# threshold<-c(1:24)
	nameFile<-paste(nameMeasure,".png",sep="")
	# tiff(file = nameFile, width = 3840, height = 2160, units = "px", res = 300)
	plot(threshold,matrix[1,], type="l", col="blue",xlab="Threshold",ylab=nameMeasure,ylim=c(0,matrix[which.max(matrix)]))
	title(main="Two-dimensional graph", col.main="red", font.main=4)
	for(i in 2:length(matrix[,1])){
		if(i<=numSeqMRNA){
			lines(threshold, matrix[i,] , type="l", pch=22, lty=2, col="blue")	
		}else{
			if(i<=(numSeqLNCRNA+numSeqMRNA)){
				lines(threshold, matrix[i,] , type="l", pch=22, lty=2, col="red")
			}else{
				lines(threshold, matrix[i,] , type="l", pch=22, lty=2, col="green")
			}	
		}
	}
	# legend((length(matrix[1,])-40), matrix[which.max(matrix)], c(1:2), cex=0.8, col=c("blue","red"), pch=21:22, lty=1:2)
	# dev.off()
	return()
}