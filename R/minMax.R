#'@title Minimum and maximum
#'@name minMax
#'
#'@description Verifies the minimum and maximum values of the results.
#'
#'@param matrix Array with results numerics
#'@param mRNA Integer number of mRNA sequences
#'@param lncRNA Integer number of lncRNA sequences
#'@param sncRNA Integer number of sncRNA sequences
#'@param rangeMinMax Vector that will be returned with the minimum and maximum values
#'
#'
#'@return Returns the vector with the minimum and maximum values for the scale
#'
#' @author Eric Augusto Ito
#'
#'

minMax <- function(matrix, mRNA, lncRNA, sncRNA, rangeMinMax){
	maxMin<-range(matrix[],na.rm = TRUE)
	rangeMinMax[2]<-maxMin[2]
	rangeMinMax[1]<-maxMin[1]
	
	if(lncRNA!=0){
		maxMin<-range(matrix[(mRNA+1):(mRNA+lncRNA),],na.rm = TRUE)
		rangeMinMax[4]<-maxMin[2]
		rangeMinMax[3]<-maxMin[1]
	}
	if(sncRNA!=0){
		maxMin<-range(matrix[(mRNA+lncRNA+1):(mRNA+lncRNA+sncRNA),],na.rm = TRUE)
		rangeMinMax[6]<-maxMin[2]
		rangeMinMax[5]<-maxMin[1]
	}
	
	return(rangeMinMax)
}