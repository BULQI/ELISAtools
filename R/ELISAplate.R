###code for defining elisa plate

#'@import minpack.lm

#defining elisa plate object, S4 class
#'@title S4 class definition of elisa_plate
#'
#' @description \code{\link{elisa_plate}} define the S4 class of elisa_plate
#'
#' @details defining the S4 class of the elisa_plate.
#'	This is the data structure to hold the elisa_plate Data. 
#' 	It contains different slots for holding both standard and 
#'    unknown data.It also defines 
#'    the regression model and the correction parameter 
#'	  for the batch effects.\cr  
#'	Note: we assume each plate has its own standard curve.
#' @concept ELISA
#'
#' @param batchID character specify batches
#'         
#' @param expID character specify experiment or plate ID
#' @param desc character data/experiment information 
#'          
#' @param data.std data.frame data for standard curves 
#'          
#' @param data.unknown data.frame data for samples with unknown concentration
##' @param model.regression nls.lm \code{link{nls.lm}} the regression model
##'		fitted with either four- or five-parameter logistic function.  
#' @param normFactor numeric the correction factor for batch effects.
#' @param mdata.unknown data.frame data to contain the mean ODs and concentration by sample IDs.
#' @param mdata.std data.frame data to contain the mean ODs and concentrations for standard data
#' @slot batchID character
#' @slot expID character
#' @slot data.std data.frame
#' @slot data.unknown data.frame
####' @slot model.regression nls.lm 
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
#'@title constructor function to build elisa_plate object
#' @description S3 method as a constructor to build the S4 
#'	class object of elisa_plate \code{\link{elisa_plate}} 
#'
#' @details S3 method as a constructor to build the S4 
#'	class object of elisa_plate \code{\link{elisa_plate}}. 
#'	Normally this is called to build a empty object 
#'	with default values and then load the data into it
#'	by calling loadData \code{\link{loadData}} or load.ODs 
#'	\code{\link{load.ODs}}
#'
#' @param batchID character specify batches
#'         
#' @param expID character specify experiment or plate ID
#' @param desc character data/experiment information 
#'          
#' @param data.std data.frame data for standard curves 
#'          
#' @param data.unknown data.frame data for samples with unknown concentration
# #' @param model.regression nls.lm \code{\link{nls.lm}} the regression model
#'		fitted with either four- or five-parameter logistic function.  
#' @param batch.normFactor numeric the correction factor for batch effects.
#' @param mdata.unknown data.frame data to contain the mean ODs and concentration by sample IDs.
#'	@return an elisa_plate object 
#' @seealso \code{\link{nls.lm}} \code{\link{loadData}} \code{\link{elisa_plate}} 
#'	\code{\link{load.ODs}}
#' @examples
#'	elisa_plate();
#' @export
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
#'@title S4 generic function to load OD data
#'
#'@description generic function to load OD data into elisa_plate 
#'		object
#'@details It load OD data into elisa_plate object. The data
#'	usually read int from a design file, annotation file and 
#'	standard concentration data.
#'
#'@param x the elisa_plate object to load data into
#'@param plate.header character
#'@param plate.data data.frame OD readings
#'@param plate.blank data.frame OD blank readings
#'@param annotation data.frame annotation to guide reading.
#'@examples
#'	elisa_plate();
#' @export 

setGeneric("load.ODs", signature="x",
			function(x, plate.header,plate.data, plate.blank,annotation,....) standardGeneric("load.ODs"))
#'@describeIn load.ODs to load ODs to the elisa_plate object
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
			ind<-paste0(rowSym[annotation$unknowns$row],annotation$unknown$col);
			x@data.unknown<-cbind(x@data.unknown, "OD"=plate.data[ind]);
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
#'@title S4 class definition of elisa_run
#'
#' @description \code{\link{elisa_run-class}} define the S4 class of elisa_run
#'
#' @details defining the S4 class of the elisa_run.
#'	This is list to hold the data for each elisa run. 
#' 	It contains one or multiple elisa plate object.  
#'
#' @concept ELISA
#'
#' @param batchID character specify batches
#'         
# # ' @param expID character specify experiment or plate ID
#' @param desc character data/experiment information 
#'          
#' @param plates list of elisa_plate 
#' @param num.plates numeric the number of plates in this run
#'
#' @slot batchID character
# #' @slot expID character
#' @slot plates list
# #' @slot data.unknown data.frame
# #' @slot model.regression nls.lm 
# #' @slot batch.normFactor numeric
#' @slot desc character
#' @slot num.plates numeric
#' @seealso \code{\link{nls.lm}} 
#' @examples
#' elisa_run();
#' @export
setClass("elisa_run",
	representation(batchID="character",#expID="character",
	desc="character", plates="list",
	num.plates="numeric", range.ODs="numeric"	
	), 
	prototype(batchID=NA_character_,#expID=NA_character_,
	desc=NA_character_, plates=list(),
	num.plates=1, range.ODs=c(-1,-1)
	)
)

#constructor
#'@title constructor function to build elisa_run object
#' @description S3 method as a constructor to build the S4 
#'	class object of elisa_run \code{\link{elisa_run-class}} 
#'
#' @details S3 method as a constructor to build the S4 
#'	class object of elisa_run \code{\link{elisa_run-class}}. 
#'	Normally this is called to build a empty object 
#'	with default values and then load the elisa_plate data into it
#'	#
#'
#' @param batchID character specify batches
#'         
#' @param desc character data/experiment information 
#'          
#' @param plates list elisa_plates in this run. could be one or multiple  
#' @param num.plates numeric the number of plates in this run.
#'
#'	@return an elisa_run object 
#' @seealso \code{\link{nls.lm}} \code{\link{elisa_run-class}} \code{\link{elisa_plate-class}} 
#'	
#' @examples
#'	elisa_run();
#' @export
elisa_run<-function(batchID=NA_character_,#expID=NA_character_,
	desc=NA_character_, plates=list(),
	num.plates=1, range.ODs=c(-1,-1)	)
{
	return(new("elisa_run",batchID=batchID,#expID=expID,
	desc=desc, plates=plates,
	num.plates=num.plates, range.ODs=range.ODs))	
}

#'@title S4 class definition of elisa_batch
#'
#' @description \code{\link{elisa_batch-class}} define the S4 class of elisa_batch
#'
#' @details defining the S4 class of the elisa_batch.
#'	This is list to hold the data for elisa batch. 
#' 	It contains one or multiple elisa_run object.  
#'
#' @concept ELISA
#'
#' @param batchID character specify batches
#'         
# # ' @param expID character specify experiment or plate ID
# @param desc character data/experiment information 
#'          
#' @param runs list of elisa_run 
#' @param num.runs numeric the number of plates in this run
#' @param pars numeric the actually parameters for the fitting.
#'		for example for 5pl it is c(a, d,xmid, scal, g).
#' @param model.fit list intend to contain information for the 
#'		fitting of nls.lm. But not using it now.
#' @param model.name character either 5pl (5-parameter) or 4pl 
#'		(4-parameter) logistic function
#' @param normFactor numeric batch-wise normalization factor 
#'
#' @slot batchID character
# #' @slot expID character
#' @slot runs list
# #' @slot data.unknown data.frame
# #' @slot model.regression nls.lm 
# #' @slot batch.normFactor numeric
#' @slot desc character
#' @slot num.runs numeric
#' @slot pars numeric 
#' @slot model.fit list
#' @slot model.name character
#' @slot normFactor numeric 
#' @seealso \code{\link{nls.lm}} \code{\link{elisa_plate-class}} \code{\link{elisa_run-class}} 
#'	\code{\link{elisa_batch-class}} 
#' @examples
#' elisa_batch();
#' @export
setClass("elisa_batch",
	representation(batchID="character",#expID="character",
	desc="character", runs="list",model.fit="list",model.name="character", pars="numeric",
	num.runs="numeric", range.ODs="numeric"	,normFactor="numeric"
	), 
	prototype(batchID=NA_character_,#expID=NA_character_,
	desc=NA_character_, runs=list(), model.fit=list(),
	model.name=NA_character_,pars=c(-1),
	num.runs=1, range.ODs=c(-1,-1), normFactor=NaN
	)
)

#constructor
#'@title constructor function to build elisa_batch object
#' @description S3 method as a constructor to build the S4 
#'	class object of elisa_batch \code{\link{elisa_batch-class}} 
#'
#' @details S3 method as a constructor to build the S4 
#'	class object of elisa_batch \code{\link{elisa_batch-class}}. 
#'	Normally this is called to build a empty object 
#'	with default values and then load the elisa_run data into it
#'	#
#'
#' @param batchID character specify batches
#'         
#' @param desc character data/experiment information 
#'          
#' @param runs list elisa_plates in this run. could be one or multiple  
#' @param num.runs numeric the number of plates in this run.
#' @param normFactor numeric batch effect normalization factor.
#'	@return an elisa_run object 
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
#'@describeIn predict to load ODs to the elisa_plate object
#'@export
setMethod("predict", c("object"="elisa_batch"),
		function(object)
		{##based on annotation, load the plate od into object
				cat("feng");
		}
)

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

#'@title predict the concentration of samples based on regression
#'@description based on the regression, predict the concentration of 
#'		of unknown sample. Assume the regression has been accomplished.
#'@details The input data structure contains both the data (ODs) and
#'		the fitted regression model. The estimation of unknonw concentration
#'		based on the OD and standard curve of each plate. And then 
#'		the batch effects are corrected/normalized and the corrected
#'		concentrations also are	also written into the batch data
#'		structure.
#'@param batches list of elisa_batch objects containing 
#'		both the raw data and the fitted regression model.
#'
#'@return the same list of elisa_batch with estimated
#'		sample concentrations based on OD and fitted regression
#'		model. The estimated concentrations normalized/corrected
#'		between different batches are also calculated and recorded.
#'@seealso  \code{\link{elisa_batch}} \code{\link{elisa_run}}
#'			\code{\link{elisa_plate}} 
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
#'@title S3 method to combine elisa_batch data
#'@description to combine the two lists of elisa_batch data.
#'     
#'@details When combining, we not only concatenate the two data sets, 
#'		but also combine batches, meaning the two 
#'		batches with same batch ID will be merge into one. 
#'		We will not merge the runs, but same batch from different
#'		list will always have different runs. It is the user's 
#'		responsibility to make sure the runs are different.  
#'@param  eb1 list of elisa_batch data 
#'@param  eb2 list of elisa_batch data 
#'
#'@return a list of elisa_batch data combining the two input lists (sorted);
#'@examples
#' setwd(system.file("extdata", package="ELISAtools"))
#' saveDB("elisa.rds");
# #' @seealso  \code{\link{elisa_batch-class}} \code{\link{loadData}} \code{\link{saveDB}}
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
				idx1<-idex1+1
			} else { #the case where eb1.ids[idx]<eb1.ids[idx2]
				batch[[count]]<-eb2[[ eb2.ids[idx2] ]];
				batch[[count]]<-resetElisaBatchAnalysis(batch[[count]]);
				idx2<-idex2+1
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
			b@runs[[i]]@plates[[j]]@data.unknown<-b@runs[[i]]@plates[[j]]@data.unknown[,c("ID","row","col","OD")]
			b@runs[[i]]@plates[[j]]@mdata.unknown<-data.frame();#b@runs[[i]]@plates[[j]]@data.std[,c("ID","conc","OD")]
		}
	}
	return(b)
}