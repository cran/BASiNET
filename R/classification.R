#'@title Performs the classification methodology using complex network theory
#'@name classification
#'
#'@description Given two distinct data sets, one of mRNA and one of lncRNA. 
#'The classification of the data is done from the structure of the networks formed by the sequences. 
#'After this is done classifying with the J48 classifier and randomForest. 
#'It is also created in the current directory a file of type arff called' result 'with all values so that it can be used later. 
#'There is also the graphic parameter that when TRUE generates graphs based on the results of each measure.
#'Using the J48 classifier it is possible to generate a tree based on the dataset and then save this tree so that it can be used to predict other RNA sequences
#'
#'@param mRNA Directory where the file .FASTA lies with the mRNA sequences
#'@param lncRNA Directory where the file .FASTA lies with the lncRNA sequences
#'@param word Integer that defines the size of the word to parse. By default the word parameter is set to 3
#'@param step Integer that determines the distance that will be traversed in the sequences for creating a new connection. By default the step parameter is set to 1
#'@param sncRNA Directory where the file .FASTA lies with the sncRNA sequences (OPTIONAL)
#'@param predicting Directory of an FASTA file containing RNA sequences. These sequences will be predicted based on the .dat file selected by the parameter "load"
#'@param graphic Parameter of the logical type, TRUE or FALSE for graphics generation. As default graphic gets FALSE
#'@param classifier Character Parameter. By default the classifier is J48, but the user can choose to use randomForest by configuring as classifier = "RF". The prediction with a model passed by the param load only works with the classifier J48.
#'@param save when set, this parameter saves a .arff file with the results of the features in the current directory and also saves the tree created by the J48 classifier so that it can be used to predict RNA sequences. This parameter sets the file name. No file is created by default
#'@param load When defined this parameter will be loaded the file which is the model previously saved in the current directory with the name entered in this parameter. No file is loaded by default
#'
#'
#' @return Results with cross-validation or the prediction result
#'
#' @author Eric Augusto Ito
#'
#'
#' @examples
#'\donttest{
#'  # Classification - cross validation
#'  library(BASiNET)
#'  arqSeqMRNA <- system.file("extdata", "sequences2.fasta", package = "BASiNET")
#'  arqSeqLNCRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
#'  classification(mRNA=arqSeqMRNA,lncRNA=arqSeqLNCRNA)
#'  classification(mRNA=arqSeqMRNA,lncRNA=arqSeqLNCRNA, save="example") #Save Tree to Predict Sequences
#'  # Prediction
#'  dataPredict <- system.file("extdata", "predict.fasta", package = "BASiNET")
#'  modelPredict <- system.file("extdata", "modelPredict.dat", package = "BASiNET")
#'  classification(predicting=dataPredict,load=modelPredict)
#'}
#' @importFrom Biostrings readBStringSet
#' @importFrom stats predict
#' @import igraph
#' @import RWeka
#' @import randomForest
#' @import rJava
#' @export

classification <- function(mRNA, lncRNA, word=3, step=1, sncRNA, predicting, graphic, classifier = c('J48', 'RF'), load, save){

	classifier <- match.arg(classifier)

	if(!missing(predicting)){
		seqMRNA<-readBStringSet(predicting)
		seqLNCRNA<-NULL
		seqSNCRNA<-NULL
		numClass<-1
	}else{
		seqMRNA<-readBStringSet(mRNA)
		seqLNCRNA<-readBStringSet(lncRNA)
		if(!missing(sncRNA)){
			seqSNCRNA<-readBStringSet(sncRNA)
			numClass<-3
		}else{
			seqSNCRNA<-NULL
			numClass<-2
		}
	}

	numSeq<-(length(seqMRNA)+length(seqLNCRNA)+length(seqSNCRNA))
	averageShortestPathLengths <- matrix(nrow=numSeq,ncol=200)
	clusteringCoefficient <- matrix(nrow=numSeq,ncol=200)
	standardDeviation <- matrix(nrow=numSeq,ncol=200)
	maximum <- matrix(nrow=numSeq,ncol=200)
	assortativity<- matrix(nrow=numSeq,ncol=200)
	betweenness <- matrix(nrow=numSeq,ncol=200)
	degree <- matrix(nrow=numSeq,ncol=200)
	minimum <- matrix(nrow=numSeq,ncol=200)
	motifs3 <- matrix(nrow=numSeq,ncol=200)
	motifs4 <- matrix(nrow=numSeq,ncol=200)

	for(k in seq_len(numClass)){
		
		if(k==1){
			if(!missing(predicting)){
				message("Analyzing RNA from number: ")
			}else{
				message("Analyzing mRNA from number: ")
			}
			aux<-0
			seq<-seqMRNA
		}else{
			if(k==2){
				message("Analyzing lncRNA from number: ")
				seq<-seqLNCRNA
				aux<-length(seqMRNA)
			}else{
				if(k==3){
					message("Analyzing sncRNA from number: ")
					seq<-seqSNCRNA
					aux<-(length(seqMRNA)+length(seqLNCRNA))
				}
			}
		}

		for(x in seq_along(seq)){
			message(x)
			sequence<-strsplit(toString(seq[x]),split='')
			sequence<-sequence[[1]]
			net<-createNet(word, step, sequence)
			limitThreshold<-max(net[])
			if(limitThreshold>200){
				limitThreshold<-200;
			}

			vector <- sapply(seq_len(limitThreshold), function(t) {
				net<<-threshold(t, net)
				measures(net)
			})
			
			cidx <- seq_len(ncol(vector))
			averageShortestPathLengths[aux + x, cidx] <- vector[1,]
			clusteringCoefficient[aux+x, cidx] <- vector[2,]
			degree[aux+x,cidx] <- vector[3,] 
			assortativity[aux+x,cidx] <- vector[4,] 
			betweenness[aux+x,cidx] <- vector[5,] 
			standardDeviation[aux+x,cidx] <- vector[6,] 
			maximum[aux+x,cidx] <- vector[7,] 
			minimum[aux+x,cidx] <- vector[8,] 
			motifs3[aux+x,cidx] <- vector[9,] 
			motifs4[aux+x,cidx] <- vector[10,] 
		}
	}

	if(!missing(load)){
		load(load)
	}

	rangeMinMax <- matrix(nrow=10,ncol=6)
	rangeMinMax[1,]<-minMax(averageShortestPathLengths,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[1,])
	rangeMinMax[2,]<-minMax(clusteringCoefficient,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[2,])
	rangeMinMax[3,]<-minMax(standardDeviation,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[3,])
	rangeMinMax[4,]<-minMax(maximum,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[4,])
	rangeMinMax[5,]<-minMax(assortativity,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[5,])
	rangeMinMax[6,]<-minMax(betweenness,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[6,])
	rangeMinMax[7,]<-minMax(degree,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[7,])
	rangeMinMax[8,]<-minMax(minimum,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[8,])
	rangeMinMax[9,]<-minMax(motifs3,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[9,])
	rangeMinMax[10,]<-minMax(motifs4,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA), rangeMinMax[10,])
	rangeMinMax[is.na(rangeMinMax)]<-0

	message("Rescaling values")
	numCol<-length(averageShortestPathLengths[1,])
	averageShortestPathLengths<-reschedule(averageShortestPathLengths,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[1,])
	colnames(averageShortestPathLengths) <- paste('ASPL',seq_len(numCol))
	clusteringCoefficient<-reschedule(clusteringCoefficient,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[2,])
	colnames(clusteringCoefficient) <- paste('CC', seq_len(numCol))
	standardDeviation<-reschedule(standardDeviation,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[3,])
	colnames(standardDeviation) <- paste('SD', seq_len(numCol))
	maximum<-reschedule(maximum,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[4,])
	colnames(maximum) <- paste('MAX', seq_len(numCol))
	assortativity<-reschedule(assortativity,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[5,])
	colnames(assortativity) <- paste('ASS', seq_len(numCol))
	betweenness<-reschedule(betweenness,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[6,])
	colnames(betweenness) <- paste('BET', seq_len(numCol))
	degree<-reschedule(degree,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[7,])
	colnames(degree) <- paste('DEG', seq_len(numCol))
	minimum<-reschedule(minimum,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[8,])
	colnames(minimum) <- paste('MIN', seq_len(numCol))
	motifs3<-reschedule(motifs3,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[9,])
	colnames(motifs3) <- paste('MT3', seq_len(numCol))
	motifs4<-reschedule(motifs4,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA),rangeMinMax[10,])
	colnames(motifs4) <- paste('MT4', seq_len(numCol))

	listMatrix<-list(
	    averageShortestPathLengths,clusteringCoefficient,
	    standardDeviation,maximum,assortativity,betweenness,degree,
	    minimum,motifs3,motifs4)
	namesMeasure<-c(
	    "Average shortest path length", "Cluster Coefficient",
	    "Standard deviation", "Maximum", "Assortativity",
	    "Betweenness", "Degree", "Minimal", "Motifs 3", "Motifs 4"
	)

	if(missing(graphic)||graphic==FALSE){
	}else{
		if(graphic==TRUE){
			for(i in seq_len(10)){
				createGraph2D(listMatrix[[i]],length(seqMRNA),length(seqLNCRNA), namesMeasure[i])
			}
		}
	}

	message("Creating data frame")
	data<-cbind(assortativity, betweenness, averageShortestPathLengths, clusteringCoefficient, degree, minimum, maximum, standardDeviation, motifs3, motifs4)
	data<-data.frame(data)
	if(numClass==2){
		data["CLASS"]<-factor(c("lncRNA"), levels = c("mRNA","lncRNA"));
		for(i in seq_along(seqMRNA)){
			data$CLASS[i] <- "mRNA"
		}
	}else{
		if(numClass==3){
			data["CLASS"]<-factor(c("lncRNA"), levels = c("mRNA","lncRNA","sncRNA"));
			data$CLASS[seq_along(seqMRNA)] <- "mRNA"
			data$CLASS[(length(seqMRNA)+length(seqLNCRNA)+1):(length(seqMRNA)+length(seqLNCRNA)+length(seqSNCRNA))] <- "sncRNA"
		}
	}
	data[is.na(data)] <- 0
	if(!missing(save)){
		rmcfs::write.arff(data, file = "Result.arff")
		message("Result.arff file generated in the current R directory")
	}

	if(classifier=="J48"){
		if(missing(load)){
			message("Sorting the data with the J48")
			obj <- J48(CLASS ~ ., data = data)
			if(!missing(save)){
				.jcache(obj$classifier)
				save(obj,file=paste(save,".dat",sep=""))
			}
			result <- evaluate_Weka_classifier(obj, numFolds = 10, complexity = TRUE, seed = 1, class = TRUE)
			print(obj)
			print(result) 
		}else{
			predict_res <- predict(object = obj, newdata=data) 
			print(predict_res)
		}
	}

	if(classifier=="RF"){
		message("Sorting the data with the Random Forest")
	    set.seed(1)
		rf <- randomForest(data[,seq_along(data[1,])], data[,"CLASS"])
		print(rf)
		print(getTree(randomForest(data[,-40], data[,5], ntree=10), 3, labelVar=TRUE))
	}

	if(missing(load)){
		return(invisible(result))
	}else{
		return(invisible(predict_res))
	}
}