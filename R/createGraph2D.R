#'@title Creates a two-dimensional graph between a measure and the threshold
#'@name createGraph2D
#'
#'@description For an analysis of each measure, the createGraph2D () function was created in order to visualize the behavior of each measurement in relation to the threshold. This function creates a graph (Measure x Threshold) from an array, mRNA sequences are given the blue color, the lncRNA sequences are given a red color. In cases where there is a third class this will be given the green color
#'
#'@param matrix matrix of the measure for the creation of two-dimensional graph
#'@param numSeqMRNA Integer number of mRNA sequences
#'@param numSeqLNCRNA Integer number of lncRNA sequences
#'@param nameMeasure Character Parameter that defines the name of the measure to put in the title of the graph
#'
#'@author Eric Augusto Ito
#'
#'
#'@import igraph
#'@importFrom grDevices dev.new
#'@importFrom graphics lines plot title


createGraph2D <- function(matrix, numSeqMRNA,numSeqLNCRNA, nameMeasure){

	dev.new()
	somador<-(1/(length(matrix[1,])-1))
	threshold<-seq(0,1,somador)
	nameFile<-paste(nameMeasure,".png",sep="")
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

}