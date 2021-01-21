###code for defining elisa plate

#'@import minpack.lm
#'@import stats
#'@importFrom methods new

#defining elisa plate object, S4 class
#'@title S4 class definition of an elisa_plate object
#'
#' @description \code{\link{elisa_plate}} define the S4 class of an elisa_plate object
#'
#' @details defining the S4 class of the elisa_plate object.
#'	This is the data structure to hold the elisa_plate Data. 
#' 	It contains different slots for holding both standard and 
#'    unknown data. It also defines 
#'    the regression model and the correction parameter 
#'	  for the batch effects.\cr  
#'	Note: we assume each plate has its own standard curve.
#' @concept ELISA
#'
#' @param batchID characters to specify the batch
#'         
#' @param expID characters to specify experiment or plate ID
#' @param desc characters for the data/experiment information 
#'          
#' @param data.std data.frame for standard curve data 
#'          
#' @param data.unknown data.frame containing data for samples with unknown concentration
# #' @param model.regression nls.lm \code{link{nls.lm}} the regression model
# #'		fitted with either the four- or five-parameter logistic function.  
#' @param normFactor numeric the correction factor for batch effects ("S").
#' @param mdata.unknown data.frame containing the mean ODs and concentration by sample IDs.
#' @param mdata.std data.frame containing the mean ODs and concentrations for standard data
#=================
#' @slot batchID character
#' @slot expID character
#' @slot data.std data.frame
#' @slot data.unknown data.frame
# ##' @slot model.regression nls.lm 
#' @slot normFactor numeric
#' @slot desc character
#' @slot range.ODs numeric
#' @slot mdata.unknown data.frame 
#' @slot mdata.std data.frame
# #' @slot offset numeric
#' @seealso \code{\link{nls.lm}} 
#' @examples
#' elisa_plate();
#' @export
setClass("elisa_plate",
	representation(batchID="character",expID="character",
	desc="character", data.std="data.frame",
	data.unknown="data.frame", #model.regression="list",
	mdata.unknown="data.frame", mdata.std="data.frame",
	normFactor="numeric", range.ODs="numeric"	
	), 
	prototype(batchID=NA_character_,expID=NA_character_,
	desc=NA_character_, data.std=data.frame(),mdata.std=data.frame(),
	data.unknown=data.frame(), mdata.unknown=data.frame(),#model.regression=list(),
	normFactor=NaN, range.ODs=c(-1,-1)	
	)
)
#constructor
#'@title Constructor function to build an elisa_plate object
#' @description S3 method as a constructor to build the S4 
#'	class object of the elisa_plate \code{\link{elisa_plate}} 
#'
#' @details S3 method as a constructor to build the S4 
#'	class object of elisa_plate \code{\link{elisa_plate}}. 
#'	Normally this is called to build an empty object 
#'	with default values and then load data into it
#'	by calling loadData \code{\link{loadData}} or load.ODs 
#'	\code{\link{load.ODs}}
#'
#' @param batchID characters to specify the batch
#'         
#' @param expID characters to specify experiment or plate ID
#' @param desc characters for the data/experiment information 
#'          
#' @param data.std data.frame for standard curve data 
#'          
#' @param data.unknown data.frame for data of samples with unknown concentration
# #' @param model.regression nls.lm \code{\link{nls.lm}} the regression model
#'		fitted with either four- or five-parameter logistic function.  
#' @param normFactor numeric the correction factor for batch effects.
#' @param mdata.unknown data.frame containing the mean ODs and concentration by sample IDs.
#' @param mdata.std data.frame containing the mean ODs and concentration of the calibration data
#' @param range.ODs numeric the min and max ODs in the plate.
#'	@return an elisa_plate object 
#' @seealso \code{\link{nls.lm}} \code{\link{loadData}} \code{\link{elisa_plate}} 
#'	\code{\link{load.ODs}}
#' @examples
#'	elisa_plate();
#' @export
# #' @importFrom methods new
elisa_plate<-function(batchID=NA_character_,expID=NA_character_,
	desc=NA_character_, data.std=data.frame(), mdata.std=data.frame(),
	data.unknown=data.frame(), mdata.unknown=data.frame(),#model.regression=list(),
	normFactor=NaN, range.ODs=c(-1,-1)	)
{
	return(new("elisa_plate",batchID=batchID,expID=expID,
	desc=desc, data.std=data.std, mdata.std=mdata.std,
	data.unknown=data.unknown, #model.regression=model.regression,
	mdata.unknown=mdata.unknown,
	normFactor=normFactor,range.ODs=range.ODs))	
}
#'@title Function to load OD data
#'
#'@description Generic function to load OD data into an elisa_plate 
#'		object
#'@details It loads OD data into an elisa_plate object. The data
#'	usually read int from design file, annotation file, OD file and 
#'	standard concentration data.
#'
#'@param x the elisa_plate object to load data into
#'@param plate.header characters
#'@param plate.data data.frame OD readings
#'@param plate.blank data.frame OD blank readings
#'@param annotation data.frame annotation to guide reading.
#'@param .... other parameters that will help reading data.
# # '@examples
# # '	x<-elisa_plate();
#' @export 

setGeneric("load.ODs", signature="x",
			function(x, plate.header,plate.data, plate.blank,annotation,....) standardGeneric("load.ODs"))
#'@describeIn load.ODs to load ODs to an elisa_plate object
setMethod("load.ODs", c("x"="elisa_plate"),
		function(x,plate.header,plate.data, plate.blank,annotation)
		{##based on annotation, load the plate od into object
			#check the data integrity 
			if(missing(x))
			{
				x<-elisa_plate();
			}
			if(missing(plate.header)||class(plate.header)!="character")
			{
				stop("please the input plate header in correct format");
			}
			if(missing(plate.data)||class(plate.data)!="numeric")
			{
				stop("please the input plate data in correct format");
			}
			if(missing(annotation)||class(annotation)!="list")
			{
				stop("please the input plate annotation in correct format");
			}
			if(missing(plate.blank)||class(plate.blank)!="numeric")
			{
				stop("please the input plate blank in correct format");
			}
			#cat("------in side load od");
			
			#now let's finally load data
			names(plate.data)<-plate.header;
			#cat("\t ^^^plate.header length:",length(plate.header),";\n")
			#cat("\t ^^^plate.blank length:",length(plate.blank),";\n")
			#cat("\t",plate.blank,"\n")
			names(plate.blank)<-plate.header;
			#cat("******done\n")
			plate.data<-plate.data-plate.blank;
			#now rewrite the row and col of annotation
			x@data.std<-annotation$standards
			x@data.unknown<-annotation$unknowns
			rowSym<-c("A","B","C","D","E","F","G","H");
			ind<-paste0(rowSym[annotation$standards$row],annotation$standards$col)
			x@data.std<-cbind(x@data.std,"OD"=plate.data[ind])
			if(!is.null(annotation$unknowns)&&dim(annotation$unknowns)[1]!=0){
				ind<-paste0(rowSym[annotation$unknowns$row],annotation$unknown$col);
				x@data.unknown<-cbind(x@data.unknown, "OD"=plate.data[ind]);
			}
			x@data.std$conc <-as.numeric(x@data.std$conc)
			#cat("====done");
			return(x)
		}#end of function
	)
	
####we also in here define two other classes of object.
# elisa_run
# elisa_batch.
#they both lists. We wraper up them so to easily work them out.
#
#defining elisa plate object, S4 class
#'@title S4 class definition of an elisa_run object
#'
#' @description \code{\link{elisa_run-class}} defines the S4 class of elisa_run
#'
#' @details defining the S4 class of the elisa_run object.
#'	This is list to hold the data for each elisa run. 
#' 	It contains one or many elisa plate objects.  
#'
#' @concept ELISA
#'
#' @param batchID characters to specify the batch
#'         
# # ' @param expID character specify experiment or plate ID
#' @param desc characters for the data/experiment information 
#'          
#' @param plates list of elisa_plates 
#' @param num.plates numeric the number of plates in this run
#'@param date characters for the date of running the ELISA measurements
#' @slot batchID character
# #' @slot expID character
#' @slot plates list
# #' @slot data.unknown data.frame
# #' @slot model.regression nls.lm 
# #' @slot batch.normFactor numeric
#' @slot desc character
#' @slot num.plates numeric
#' @slot date character
#' @slot range.ODs numeric 
#' @seealso \code{\link{nls.lm}} 
#' @examples
#' elisa_run();
#' @export
setClass("elisa_run",
	representation(batchID="character",#expID="character",
	desc="character", plates="list",
	num.plates="numeric", date="character",range.ODs="numeric"	
	), 
	prototype(batchID=NA_character_,#expID=NA_character_,
	desc=NA_character_, plates=list(),
	num.plates=0, date="NA_character_",range.ODs=c(-1,-1)
	)
)

#constructor
#'@title Constructor function to build an elisa_run object
#' @description S3 method as a constructor to build the S4 
#'	class object of elisa_run \code{\link{elisa_run-class}} 
#'
#' @details S3 method as a constructor to build the S4 
#'	class object of elisa_run \code{\link{elisa_run-class}}. 
#'	Normally this is called to build a empty object 
#'	with default values and then load the elisa_run data into.
#'	#
#'
#' @param batchID characters to specify the batch
#'         
#' @param desc characters for data/experiment information 
#'          
#' @param plates list of elisa_plates in this run. could be one or many  
#' @param num.plates numeric the number of plates in this run.
#' @param date charaters the date to run ELISA measurements
#' @param range.ODs numeric the range of ODs for the measurements
#'	@return an elisa_run object 
#' @seealso \code{\link{nls.lm}} \code{\link{elisa_run-class}} \code{\link{elisa_plate-class}} 
#'	
#' @examples
#'	elisa_run();
#' @export
elisa_run<-function(batchID=NA_character_,#expID=NA_character_,
	desc=NA_character_, plates=list(),
	num.plates=1, date=NA_character_, range.ODs=c(-1,-1)	)
{
	return(new("elisa_run",batchID=batchID,#expID=expID,
	desc=desc, plates=plates,
	num.plates=num.plates, date=date,range.ODs=range.ODs))	
}

#'@title S4 class definition of an elisa_batch object
#'
#' @description \code{\link{elisa_batch-class}} define the S4 class of an elisa_batch object
#'
#' @details defining the S4 class of the elisa_batch.
#'	This holds the data for elisa batch. 
#' 	It contains one or many elisa_run objects.  
#'
#' @concept ELISA
#'
#' @param batchID characters to specify the batch
#'         
# # ' @param expID character specify experiment or plate ID
# @param desc character string for the data/experiment information 
#'          
#' @param runs list of elisa_run objects 
#' @param num.runs numeric the number of elisa_runs in this batch
#' @param pars numeric the actually parameters for the fitting.
#'		for example for the 5pl they are c(a, d,xmid, scal, g).
#' @param model.fit list intend to contain information for the 
#'		fitting of nls.lm. But not using it now.
#' @param model.name characters of either the 5pl (5-parameter) or 4pl 
#'		(4-parameter) logistic function
#' @param range.ODs numeric the min and max ODs
#' @param normFactor numeric the batch normalization factor ("S"). 
#'
#==========================
#' @slot batchID character
# #' @slot expID character
#' @slot desc character
#' @slot runs list
#' @slot num.runs numeric
# #' @slot data.unknown data.frame
# #' @slot model.regression nls.lm 
# #' @slot batch.normFactor numeric
#' @slot pars numeric 
#' @slot model.fit list
#' @slot model.name character
#' @slot range.ODs numeic
#' @slot normFactor numeric 
#' @seealso \code{\link{nls.lm}} \code{\link{elisa_plate-class}} \code{\link{elisa_run-class}} 
#'	\code{\link{elisa_batch-class}} 
#' @examples
#' elisa_batch();
#' @export
setClass("elisa_batch",
	representation(batchID="character",#expID="character",
	desc="character", runs="list", num.runs="numeric", pars="numeric",
	model.fit="list",model.name="character", 
	range.ODs="numeric"	,normFactor="numeric"
	), 
	prototype(batchID=NA_character_,#expID=NA_character_,
	desc=NA_character_, runs=list(), model.fit=list(),
	model.name=NA_character_,pars=c(-1),
	num.runs=1, range.ODs=c(-1,-1), normFactor=NaN
	)
)

#constructor
#'@title Constructor function to build an elisa_batch object
#' @description S3 method as a constructor to build the S4 
#'	class object of elisa_batch \code{\link{elisa_batch-class}} 
#'
#' @details S3 method as a constructor to build the S4 
#'	class object of elisa_batch \code{\link{elisa_batch-class}}. 
#'	Normally this is called to build a empty object 
#'	with default values and then load the elisa_run data into it
#'
#'
#' @param batchID character string to specify a batches
#'         
#' @param desc character string for the data/experiment information 
#'          
#' @param runs list of elisa_plates in this run. There could be one or many plates in a run.  
#' @param num.runs numeric the number of plates in this run.
#' @param pars numeric the actually parameters for the fitting.
#'		for example for the 5pl model they are c(a, d,xmid, scal, g).
#' @param model.fit list intend to contain information for the 
#'		fitting of nls.lm. But not using it now.
#' @param model.name character string of either the 5pl (5-parameter) or 4pl 
#'		(4-parameter) logistic function
#' @param range.ODs numeric the min and max ODs
#'
#' @param normFactor numeric the batch effect normalization factor ("S").
#'	@return an elisa_batch object 
#' @seealso \code{\link{nls.lm}} \code{\link{elisa_run-class}} \code{\link{elisa_plate-class}} 
#'	\code{\link{elisa_batch-class}} 
#' @examples
#'	elisa_batch();
#' @export
elisa_batch<-function(batchID=NA_character_,#expID=NA_character_,
	desc=NA_character_, runs=list(), model.fit=list(),model.name=NA_character_,pars=c(-1),
	num.runs=1,	range.ODs=c(-1,-1),normFactor=NaN)
{
	return(new("elisa_batch",batchID=batchID,#expID=expID,
	desc=desc, runs=runs,model.fit=model.fit, model.name=model.name, pars=pars,
	num.runs=num.runs, range.ODs=range.ODs, normFactor=normFactor))	
}

#setClassUnion("nls.lmOrNULL",c("nls.lm","NULL")) 
#setGeneric
#setGeneric("predict", signature="x",
#			function(x, plate.header,plate.data, plate.blank,annotation,....) standardGeneric("load.ODs"))
# #'@describeIn predict to load ODs to the elisa_plate object
# #'@export
#setMethod("predict", c("object"="elisa_batch"),
#		function(object)
#		{##based on annotation, load the plate od into object
#				cat("feng");
#		}
#)

#predict based on the fitted model the elisa_batch data
#either with batch correction or without batch correction 
# #'@export
predictBatchData<-function(batch)
{
	if(missing(batch))
	{
		stop("ERROR:please specify the input batch data");
	}
	if(class(batch)!="elisa_batch")
	{
		stop("ERROR:please specify the data as elisa_batch");
	}
	
	#now get the model
	pars<-batch@pars
	batchNormFac<-batch@normFactor;
	
	#predict.mean<-list();
	#now go through each 
	count<-0;
	for(i in 1:batch@num.runs)
	{
		for(j in 1:batch@runs[[i]]@num.plates)
		{
			count<-count+1;
			#get the model
			pars.plate<-pars+c(0,0,-1*batch@runs[[i]]@plates[[j]]@normFactor,0,0)
			
			#for doing unknows
			if(!is.null(batch@runs[[i]]@plates[[j]]@data.unknown)&& dim(batch@runs[[i]]@plates[[j]]@data.unknown)[1]!=0){
				x.expect<-inv.f5pl(pars.plate, batch@runs[[i]]@plates[[j]]@data.unknown$OD);
				batch@runs[[i]]@plates[[j]]@data.unknown$conc_pred<-x.expect;
				batch@runs[[i]]@plates[[j]]@data.unknown$conc_pred.bc<-x.expect*exp(batchNormFac);
				
				#doing 
				dF<-aggregate(batch@runs[[i]]@plates[[j]]@data.unknown$OD, 
						by=list(group=batch@runs[[i]]@plates[[j]]@data.unknown$ID), FUN=mean);
				
				x.expect<-inv.f5pl(pars.plate, dF$x);
				dF$conc_pred<-x.expect
				dF$conc_pred.bc<-x.expect*exp(batchNormFac);
				dF<-dF[,c("group","x","conc_pred","conc_pred.bc")]
				names(dF)<-c("ID","OD","conc_pred","conc_pred.bc")
				batch@runs[[i]]@plates[[j]]@mdata.unknown<-dF;
			}
			#--now add the predications for standard table
			stdx.expect<-inv.f5pl(pars.plate, batch@runs[[i]]@plates[[j]]@data.std$OD);
			batch@runs[[i]]@plates[[j]]@data.std$conc_pred<-stdx.expect;
			batch@runs[[i]]@plates[[j]]@data.std$conc_pred.bc<-stdx.expect*exp(batchNormFac);
			dF<-aggregate(batch@runs[[i]]@plates[[j]]@data.std$OD, 
					by=list(group=batch@runs[[i]]@plates[[j]]@data.std$ID), FUN=mean);
			dF.conc<-aggregate(batch@runs[[i]]@plates[[j]]@data.std$conc, 
					by=list(group=batch@runs[[i]]@plates[[j]]@data.std$ID), FUN=mean);
			names(dF.conc)<-c("ID","conc");
			stdx.expect<-inv.f5pl(pars.plate, dF$x);
			dF$conc_pred<-stdx.expect
			dF$conc_pred.bc<-stdx.expect*exp(batchNormFac);
			dF<-dF[,c("group","x","conc_pred","conc_pred.bc")]
			names(dF)<-c("ID","OD","conc_pred","conc_pred.bc")
			dF.conc<-cbind(dF.conc,dF[,c("OD","conc_pred","conc_pred.bc")]);
			batch@runs[[i]]@plates[[j]]@mdata.std<-dF.conc;
			
		}
	}
	return (batch)
}

#'@title Predict the concentration of samples based on fitting
#'@description Based on the 5pl or 4pl regression, predict the concentration of 
#'		of unknown samples. Assume the regression has been accomplished.
#'@details The input data structure contains both the data (ODs) and
#'		the fitted regression model. The estimation of unknonw concentration
#'		based on the ODs and the standard curve of each plate.  
#'		The batch effects are corrected/normalized and the corrected
#'		concentrations also are	also written into the batch data
#'		structure, if there are more than one batches in the data.
#'@param batches list of elisa_batch objects containing 
#'		both the raw data and the fitted regression model.
#'
#'@return The same list of elisa_batch with estimated
#'		sample concentrations based on ODs and the fitted regression
#'		model. The estimated concentrations normalized/corrected
#'		between different batches are also calculated and recorded.
#'@seealso  \code{\link{elisa_batch}} \code{\link{elisa_run}}
#'			\code{\link{elisa_plate}} 
#'@references Feng 2018 \doi{10.1101/483800}
#'@export
predictAll<-function(batches)
{
	if(missing(batches))
	{
		stop("ERROR:please specify the input batch data");
	}
	if(class(batches)!="list")
	{
		stop("ERROR:please specify the data as a list of elisa_batch object");
	}
	
	#now go through the data to get everything done
	
	for(i in 1:length(batches))
	{
		batch<-batches[[i]]
		batches[[i]]<-predictBatchData(batch)
	}
	return(batches)
}
#functions to combine two elisa batch list into one.

#----------for data base saving/loading
#functions to 
#'@title Combine elisa_batch data
#'@description Combine the two lists of elisa_batch data.
#'     
#'@details When combining, we not only concatenate the two data sets, 
#'		but also combine batches, meaning the two 
#'		batches with same batch ID will be merged into one. 
#'		We will not merge the runs. Therefore, same batch from different
#'		list will always have different runs. It is the user's 
#'		responsibility to make sure the runs are different.  
#'@param  eb1 list of elisa_batch data 
#'@param  eb2 list of elisa_batch data 
#'
#'@return a list of elisa_batch data combining the two input lists (sorted);
#'@examples
#'#R code to run 5-parameter logistic regression on ELISA data
#'#load the library
#'library(ELISAtools)
#'
#'#get file folder
#'dir_file<-system.file("extdata", package="ELISAtools")
#'
#'batches<-loadData(file.path(dir_file,"design.txt"))
#'
#'#make a guess for the parameters, the other two parameters a and d 
#'#will be estimated based on data.
#'model<-"5pl"
#'pars<-c(7.2,0.5, 0.015) #5pl inits
#'names(pars)<-c("xmid", "scal", "g")
#'
#'#do fitting. model will be written into data set.
#'batches<-runFit(pars=pars,  batches=batches, refBatch.ID=1, model=model  )
#'
#'#call to do predications based on the model.
#'batches<-predictAll(batches);
#'
#'batches.old<-batches;
#'
#'#now suppose want to join/combine the two batches, old and new 
#'batches.com<-combineData(batches.old, batches);
#'
#'@seealso  \code{\link{elisa_batch-class}} \code{\link{loadData}} \code{\link{saveDB}}
#'
#'@export
combineData<-function(eb1, eb2)
{
		if(missing(eb1)||missing(eb2))
		{
			stop("please specify the input data")
		}
		#let's do merge sort kind of combining.
		#first get the batch ids
		eb1.id<-names(eb1);
		eb2.id<-names(eb2);
		
		#sort them
		eb1.ids<-sort(eb1.id)
		eb2.ids<-sort(eb2.id)
		
		#now go through to combine then
		len1<-length(eb1.ids)
		len2<-length(eb2.ids)
		idx1<-1;
		idx2<-1;
		count<-0;
		#flag<-TRUE
		batch<-list();
		while(TRUE){
			#check for orders
			if(idx1>len1 || idx2>len2)
			{
				#flag<-FALSE;
				break;
			}
			count<-count+1;
			if(eb1.ids[idx1]==eb2.ids[idx2])
			{
				#combine batch together.
				batch[[count]]<-combineBatch(eb1[[eb1.ids[idx1] ]],eb2[[eb2.ids[idx2] ]]);
				idx1<-idx1+1;
				idx2<-idx2+1;
				#next;
			} else if(eb1.ids[idx1]>eb1.ids[idx2]){
				batch[[count]]<-eb1[[ eb1.ids[idx1] ]];
				batch[[count]]<-resetElisaBatchAnalysis(batch[[count]]);
				idx1<-idx1+1
			} else { #the case where eb1.ids[idx]<eb1.ids[idx2]
				batch[[count]]<-eb2[[ eb2.ids[idx2] ]];
				batch[[count]]<-resetElisaBatchAnalysis(batch[[count]]);
				idx2<-idx2+1
			}
		}#end of merge combine.
		
		#now need to copy over the left-behind
		if(idx1<=len1)
		{
			for(i in idx1:len1)
			{
				count<-count+1;
				batch[[count]]<-eb1[[ eb1.ids[i] ]]
			}
		}
		if(idx2<=len2)
		{
			for(i in idx2:len2)
			{
				count<-count+1;
				batch[[count]]<-eb2[[ eb2.ids[i] ]];
			}
		}
		return(batch)
}
#internal function
#we will combine two batches, with identical batch IDs
#	When combined, we will assume the two have identical
#	batchIDs and will combine the runs, reset the analysis
#	assuming we will have to do the analysis again.
# #'@param b1 elisa_batch data
# #'@param b2 elisa_batch data
# #'@return combined elisa_batch data.

combineBatch<-function(b1,b2)
{
	batch<-b1;
	#copy over the runs
	num.runs2<-b2@num.runs;
	num.runs1<-b1@num.runs;
	for(i in 1:num.runs2)
	{
		batch@runs[[num.runs1+i]]<-b2@runs[[i]];
	}
	batch@num.runs<-num.runs1+num.runs2;
	batch@range.ODs[1]<-min(b1@range.ODs[1], b2@range.ODs[1]);
	batch@range.ODs[2]<-max(b1@range.ODs[2], b2@range.ODs[2]);
	
	#always reset the analysis.
	#if(is.na(b1@normFactor)||is.na(b2@normFactor)||b1@normFactor!=b2@normFactor)
	#{
	resetElisaBatchAnalysis(batch);
	#}
	#
	#return(batch)
}

#internal function to reset the analysis of
#the elisa_batch data.
resetElisaBatchAnalysis<-function(b)
{
	b@normFactor<- NaN;
	b@pars<-c(-1);
	b@model.name<-"";
	for(i in 1:b@num.runs)
	{
		for(j in 1:b@runs[[i]]@num.plates)
		{
			#reset the predication in both tables
			b@runs[[i]]@plates[[j]]@normFactor<-NaN
			b@runs[[i]]@plates[[j]]@data.std<-b@runs[[i]]@plates[[j]]@data.std[,c("ID","row","col","conc","OD")]
			b@runs[[i]]@plates[[j]]@mdata.std<-data.frame();#b@runs[[i]]@plates[[j]]@data.std[,c("ID","conc","OD")]
			if(!is.null(b@runs[[i]]@plates[[j]]@data.unknown)&&dim(b@runs[[i]]@plates[[j]]@data.unknown)[1]!=0){
				b@runs[[i]]@plates[[j]]@data.unknown<-b@runs[[i]]@plates[[j]]@data.unknown[,c("ID","row","col","OD")];
			}	
			b@runs[[i]]@plates[[j]]@mdata.unknown<-data.frame();#b@runs[[i]]@plates[[j]]@data.std[,c("ID","conc","OD")]
		}
	}
	return(b)
}