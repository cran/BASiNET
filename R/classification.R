#'@title Performs the classification methodology using complex network theory
#'@name classification
#'
#'@description Given two distinct data sets, one of mnRNA and one of lncRNA. 
#'The classification of the data is done from the structure of the networks formed by the sequences. 
#'After this is done classifying with the J48 classifier and randomForest. 
#'It is also created in the current directory a file of type arff called' result 'with all values so that it can be used later. 
#There is also the graphic parameter that when TRUE generates graphs based on the results of each measure.
#'
#'@param mRNA Directory where the file lies with the mRNA sequences
#'@param lncRNA Directory where the file is located fasta with lncRNA sequences
#'@param word Size of the word to parse. By default the word parameter is set to 3
#'@param step Determines the distance that will be traversed in the sequences for creating a new connection. By default the step parameter is set to 1
#'@param sncRNA Directory where the file is located fasta with the sncRNA sequences (OPTIONAL)
#'@param graphic TRUE or FALSE for graphics generation. As default graphic gets FALSE
#'@param graphic3D TRUE or FALSE for 3D graphics generation. As default graphic3D gets FALSE
#'@param classifier By default the classifier is J48, but the user can choose to use randomForest by configuring as classifier = "RF"
#'
#'
#' @return Data.frame with the results of measures
#'
#' @author Eric Augusto Ito
#'
#'
#' @examples
#' arqSeqMRNA <- system.file("extdata", "sequences2.fasta", package = "BASiNET")
#' arqSeqLNCRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
#' classification(mRNA=arqSeqMRNA,lncRNA=arqSeqLNCRNA,word=3,step=3,graphic=FALSE,graphic3D=FALSE)
#'
#' @importFrom Biostrings readBStringSet
#' @import igraph
#' @import RWeka
#' @import randomForest
#' @export

classification <- function(mRNA, lncRNA, word, step, sncRNA, graphic, graphic3D, classifier){


	if(missing(classifier)){
		classifier<-"J48"
	}

	if((classifier!="RF")&&(classifier!="J48")){
		print("INVALID CLASSIFIER! Classifier configured for J48")
		classifier<-"J48"
	}
	if(missing(word)){
		word<-3
	}
	if(missing(step)){
		step <- 1
	}

	seqMRNA<-readBStringSet(mRNA)
	seqLNCRNA<-readBStringSet(lncRNA)

	if(!missing(sncRNA)){
		seqSNCRNA<-readBStringSet(sncRNA)
		numClass<-3
	}else{
		seqSNCRNA<-NULL
		numClass<-2
	}

	numSeq<-(length(seqMRNA)+length(seqLNCRNA)+length(seqSNCRNA))
	averageShortestPathLengths <- matrix(nrow=numSeq,ncol=1)    
	clusteringCoefficient <- matrix(nrow=numSeq,ncol=1) 
	standardDeviation <- matrix(nrow=numSeq,ncol=1) 
	maximum <- matrix(nrow=numSeq,ncol=1) 
	assortativity<- matrix(nrow=numSeq,ncol=1) 
	betweenness <- matrix(nrow=numSeq,ncol=1) 
	degree <- matrix(nrow=numSeq,ncol=1) 
	minimum <- matrix(nrow=numSeq,ncol=1) 
	motifs3 <- matrix(nrow=numSeq,ncol=1) 
	motifs4 <- matrix(nrow=numSeq,ncol=1)

	for(k in 1:numClass){
		
		if(k==1){
			print("Analyzing mRNA from number: ")
			aux<-0
			seq<-seqMRNA
			type<-"mRNA"
		}else{
			if(k==2){
				print("Analyzing lncRNA from number: ")
				seq<-seqLNCRNA
				aux<-length(seqMRNA)
				type<-"lncRNA"
			}else{
				if(k==3){
					print("Analyzing sncRNA from number: ")
					seq<-seqSNCRNA
					aux<-(length(seqMRNA)+length(seqLNCRNA))
				}
			}
		}

		for(x in 1:length(seq)){
			print(x)
			limitThreshold<-0
			sequence<-strsplit(toString(seq[x]),split='')
			sequence<-sequence[[1]]
			net<-createNet(word, step, sequence)
			limitThreshold<-max(net[])

			for(t in 1:limitThreshold){
				if(t>length(averageShortestPathLengths[1,])){
					averageShortestPathLengths <- cbind(averageShortestPathLengths,matrix(nrow=numSeq,ncol=1))
					clusteringCoefficient <- cbind(clusteringCoefficient,matrix(nrow=numSeq,ncol=1))
					standardDeviation <- cbind(standardDeviation,matrix(nrow=numSeq,ncol=1))
					maximum <- cbind(maximum,matrix(nrow=numSeq,ncol=1))
					assortativity <- cbind(assortativity,matrix(nrow=numSeq,ncol=1))
					betweenness <- cbind(betweenness,matrix(nrow=numSeq,ncol=1))
					degree <-cbind(degree,matrix(nrow=numSeq,ncol=1))
					minimum <- cbind(minimum,matrix(nrow=numSeq,ncol=1))
					motifs3 <- cbind(motifs3,matrix(nrow=numSeq,ncol=1))
					motifs4 <- cbind(motifs4,matrix(nrow=numSeq,ncol=1))
				}

				net<-threshold(t, net)
				vector<-measures(net)
				averageShortestPathLengths[aux+x,t] <- vector[1]    
				clusteringCoefficient[aux+x,t] <- vector[2] 
				degree[aux+x,t] <- vector[3] 
				assortativity[aux+x,t] <- vector[4] 
				betweenness[aux+x,t] <- vector[5] 
				standardDeviation[aux+x,t] <- vector[6] 
				maximum[aux+x,t] <- vector[7] 
				minimum[aux+x,t] <- vector[8] 
				motifs3[aux+x,t] <- vector[9] 
				motifs4[aux+x,t] <- vector[10] 
				
			}
		}
	}
	print("Rescaling values")
	numCol<-length(averageShortestPathLengths[1,])
	averageShortestPathLengths<-reschedule(averageShortestPathLengths,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(averageShortestPathLengths) <- paste('ASPL',1:numCol)
	clusteringCoefficient<-reschedule(clusteringCoefficient,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(clusteringCoefficient) <- paste('CC', 1:numCol)
	standardDeviation<-reschedule(standardDeviation,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(standardDeviation) <- paste('SD', 1:numCol)
	maximum<-reschedule(maximum,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(maximum) <- paste('MAX', 1:numCol)
	assortativity<-reschedule(assortativity,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(assortativity) <- paste('ASS', 1:numCol)
	betweenness<-reschedule(betweenness,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(betweenness) <- paste('BET', 1:numCol)
	degree<-reschedule(degree,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(degree) <- paste('DEG', 1:numCol)
	minimum<-reschedule(minimum,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(minimum) <- paste('MIN', 1:numCol)
	motifs3<-reschedule(motifs3,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(motifs3) <- paste('MT3', 1:numCol)
	motifs4<-reschedule(motifs4,length(seqMRNA),length(seqLNCRNA),length(seqSNCRNA))
	colnames(motifs4) <- paste('MT4', 1:numCol)

	listMatrix<-list(averageShortestPathLengths,clusteringCoefficient,standardDeviation,maximum,assortativity,betweenness,degree,minimum,motifs3,motifs4)
	namesMeasure<-c("Average shortest path length", "Cluster Coefficient", "Standard deviation", "Maximum", "Assortativity", "Betweenness", "Degree", "Minimal", "Motifs 3", "Motifs 4")
 

	
	if(missing(graphic)||graphic==FALSE){
	}else{
		if(graphic==TRUE){
			for(i in 1:10){
				createGraph2D(listMatrix[[i]],length(seqMRNA),length(seqLNCRNA), namesMeasure[i])
			}
		}
	}	
	if(missing(graphic3D)||graphic3D==FALSE){
	}else{
		if(graphic3D==TRUE){
			for(x in 1:10){
				for(y in 1:10){
					if(y>x){
						createGraph3D(listMatrix[[x]],listMatrix[[y]],length(seqMRNA),length(seqLNCRNA),namesMeasure[x],namesMeasure[y])
					}
				}
			}
		}
	}

	
	#cria uma janela em branco para exibicao da arvore do J48
	# class1<-matrix(nrow=length(seqMRNA),ncol=1)
	# class2<-matrix(nrow=length(seqLNCRNA),ncol=1)
	# class1[,1]<-"mRNA"
	# class2[,1]<-"lncRNA"
	# classes<-rbind(class1,class2)
	# colnames(classes) <- paste('CLASS') 

	print("Creating data frame")
	data<-cbind(assortativity,betweenness,averageShortestPathLengths, clusteringCoefficient,degree,minimum,maximum,standardDeviation,motifs3,motifs4)
	data<-data.frame(data)
	if(numClass==2){
		data["CLASS"]<-factor(c("lncRNA"), levels = c("mRNA","lncRNA"));
		for(i in 1:length(seqMRNA)){
			data$CLASS[i] <- "mRNA"
		}
	}else{
		if(numClass==3){
			data["CLASS"]<-factor(c("lncRNA"), levels = c("mRNA","lncRNA","sncRNA"));
			data$CLASS[1:length(seqMRNA)] <- "mRNA"
			data$CLASS[(length(seqMRNA)+length(seqLNCRNA)+1):(length(seqMRNA)+length(seqLNCRNA)+length(seqSNCRNA))] <- "sncRNA"
		}
	}
	data[is.na(data)] <- 0
	rmcfs::write.arff(data, file = "Result.arff")
	print("Result.arff file generated in the current R directory")

	if(classifier=="J48"){
		print("Sorting the data with the J48")
		obj <- J48(CLASS ~ ., data = data)
		result <- evaluate_Weka_classifier(obj, numFolds = 10, complexity = TRUE, seed = 1, class = TRUE)
		print(obj)
		print(result)
		# plot(obj)
	}

	if(classifier=="RF"){
		print("Sorting the data with the Random Forest")
	    set.seed(1)
		rf <- randomForest(data[,1:length(data[1,])], data[,"CLASS"])
		print(rf)
		print(getTree(randomForest(data[,-40], data[,5], ntree=10), 3, labelVar=TRUE))
	}

	return(invisible(data))
}