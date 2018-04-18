#'@title Rescales the results between values from 0 to 1
#'@name reschedule
#'
#'@description Given the results the data is rescaled for values between 0 and 1, so that the length of the sequences does not influence the results. The rescaling of the mRNA and lncRNA are made separately
#'
#'@param matrix Array with results numerics
#'@param mRNA Integer number of mRNA sequences
#'@param lncRNA Integer number of lncRNA sequences
#'@param sncRNA Integer number of sncRNA sequences
#'@param rangeMinMax Vector with the minimum and maximum values for the scale
#'
#'@return Returns the array with the rescaled values
#'
#' @author Eric Augusto Ito
#'
#'

reschedule <- function(matrix, mRNA, lncRNA, sncRNA, rangeMinMax){

	for(x in 1:mRNA){
		for(y in 1:length(matrix[1,])){
			matrix[x,y]<-((matrix[x,y]-rangeMinMax[1])/(rangeMinMax[2]-rangeMinMax[1]))
		}
	}
	if(lncRNA!=0){
		for(x in (mRNA+1):(mRNA+lncRNA)){
			for(y in 1:length(matrix[1,])){
				matrix[x,y]<-((matrix[x,y]-rangeMinMax[3])/(rangeMinMax[4]-rangeMinMax[3]))
			}
		}
	}
	if(sncRNA!=0){
		for(x in (mRNA+lncRNA+1):(mRNA+lncRNA+sncRNA)){
			for(y in 1:length(matrix[1,])){
				matrix[x,y]<-((matrix[x,y]-rangeMinMax[5])/(rangeMinMax[6]-rangeMinMax[5]))
			}
		}
	}
	
	return(matrix)
}