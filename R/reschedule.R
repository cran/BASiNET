#'@title Rescales the results between values from 0 to 1
#'@name reschedule
#'
#'@description Given the results the data is rescaled for values between 0 and 1, so that the length of the sequences does not influence the results. The rescaling of the mRNA and lncRNA are made separately
#'
#'@param matrix array with results
#'@param mRNA Number of mRNA sequences
#'@param lncRNA Number of lncRNA sequences
#'@param sncRNA Number of sncRNA sequences
#'
#'
#'@return Returns the array with the rescaled values
#'
#' @author Eric Augusto Ito
#'
#'


reschedule <- function(matrix, mRNA, lncRNA, sncRNA){
	maxMin<-range(matrix[1:mRNA,],na.rm = TRUE)
	max<-maxMin[2]
	min<-maxMin[1]
	for(x in 1:mRNA){
		for(y in 1:length(matrix[1,])){
			matrix[x,y]<-((matrix[x,y]-min)/(max-min))
		}
	}

	maxMin<-range(matrix[(mRNA+1):(mRNA+lncRNA),],na.rm = TRUE)
	max<-maxMin[2]
	min<-maxMin[1]
	for(x in (mRNA+1):(mRNA+lncRNA)){
		for(y in 1:length(matrix[1,])){
			matrix[x,y]<-((matrix[x,y]-min)/(max-min))
		}
	}
	
	if(sncRNA!=0){
		maxMin<-range(matrix[(mRNA+lncRNA+1):(mRNA+lncRNA+sncRNA),],na.rm = TRUE)
		max<-maxMin[2]
		min<-maxMin[1]
		for(x in (mRNA+lncRNA+1):(mRNA+lncRNA+sncRNA)){
			for(y in 1:length(matrix[1,])){
				matrix[x,y]<-((matrix[x,y]-min)/(max-min))
			}
		}
	}
	
	return(matrix)
}