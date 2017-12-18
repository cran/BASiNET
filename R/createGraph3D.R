#'@title Creates an untargeted graph from a biological sequence
#'@name createGraph3D
#'
#'@description A function that from a biological sequence generates a graph not addressed having as words vertices, this being able to have its size parameter set by the' word 'parameter. The connections between words depend of the' step 'parameter that indicates the next connection to be formed
#'
#'@param matrix1 matrix of the first measure for the creation of the three-dimensional chart
#'@param matrix2 matrix of the second measure for the creation of the three-dimensional chart
#'@param numSeqMRNA number of mRNA sequences
#'@param numSeqLNCRNA number of lncRNA sequences
#'@param nameMeasure1 name of the first measure to put in the title of the graph
#'@param nameMeasure2 name of the second measure to put in the title of the graph
#'
#'@return Returns the non-directed graph formed through the sequence
#'
#' @author Eric Augusto Ito
#'
#' @import rgl



createGraph3D <- function(matrix1, matrix2, numSeqMRNA, numSeqLNCRNA, nameMeasure1, nameMeasure2){

	somador<-(1/(length(matrix1[1,])-1))
	threshold<-seq(0,1,somador)
	open3d()
	plot3d(matrix1[1,],matrix2[1,], threshold, xlab=nameMeasure1, ylab=nameMeasure2, zlab="Threshold", col=c("blue"), type="l2")
	for(i in 2:length(matrix1[,1])){
		if(i<=numSeqMRNA){
			lines3d(matrix1[i,], matrix2[i,], threshold,color=c("blue"))	
		}else{
			if(i<=(numSeqLNCRNA+numSeqMRNA)){
				lines3d(matrix1[i,], matrix2[i,], threshold,color=c("red"))	
			}else{
				lines3d(matrix1[i,], matrix2[i,], threshold,color=c("green"))	
			}
		}
	}

	return()
}