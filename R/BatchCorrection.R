
#'@import R2HTML
#'@import grDevices
#'@import graphics


#BatchCorretion module for ELISA analysis
#it is in this model to prepare the input and parameters
#	 and then pass them onto regression model for fitting
#We will also in this module to plotting and diagnosing, etc.
#
#in this function, we prepare input for the regression.
#we need to take into consideration of batch effects
#1) make x, y ready
#2) make initial values ready
#3) prepare a, d 
#4) prepare k, shifts., figure out the aggregation matrix, to repeat k's
#'@title Prepare the input for regressoin
#'@description  Prepare the input data to feed in the fitting.
#'
#'@param batches list of the ELISA data arranged in 
#'	batches. Each element of list contains a batch (list)
#'	data, and each batch contains one or many the elisa_run 
#'	objects \code{\link{elisa_plate}}
#'	
# #'@param ref.ID characters the reference batch ID. all other 
# #'	batches are shifted/corrected towards this reference batch
#' 
#'@return list of data that will feed in to do regression.
#'@seealso \code{\link{elisa_plate}} \code{\link{elisa_batch}}\code{\link{elisa_run}}
#'@export 
prepareRegInput<-function(batches
		)
{
	if(missing(batches))
	{
		stop("please specify the input batch data")
	}
	#go through data and get all the data into data frames, x and y, separately
	#we also record the number of points in each standard, in order to 
	#account for unequal x values. we will use this to expand 
	# k (shifts) during residual function
	ind.batch<-c();
	ind.run<-c();
	ind.plate<-c();
	num.std<-c()
	y<-c();
	x<-c()
	for(i in 1:length(batches))
	{
		ebatch<-batches[[i]]
		for(j in 1: ebatch@num.runs)
		{
			erun<-ebatch@runs[[j]]
			for(k in 1:erun@num.plates)
			{
				eplate<-erun@plates[[k]]
				ind.batch<-c(ind.batch, i)
				ind.run<-c(ind.run,j)
				ind.plate<-c(ind.plate,k)
				num.std<-c(num.std,length(eplate@data.std$OD))
				y<-c(y, eplate@data.std$OD)
				x<-c(x, eplate@data.std$conc)
			}
		}
	}
	
	#now we have everything
	return(list(y=y,x=x,ind.batch=ind.batch, ind.run=ind.run, ind.plate=ind.plate, num.stds=num.std))
}

    #'@title Prepare initial values for fitting shifts 
	#'@description Generate the initial values for fitting shifts with 
	#'	a model of the 5-parameter logistic function.
	#'@details This is a more complicated way to prepare the initials for shifting. 
	# #'	The assumption is that the OD values of analytes in ELISA 
	# #'	follow 
	# #'	a 5-parameter logistic function (5pl) as \cr
#	#'	\ifelse{html}{\out{<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML"></script>
#	#'						$$y = a + {{d-a} \over {(1 + e^ {{(xmid - x)} \over scal})}^g} $$
#	#'				 	}}{\deqn{y=a+\frac{d-a}{(1+e^(\frac{xmid-x}{scal}))^g}}{ASCII}}
#	#'	\cr 
#	#'	where a is the upper bound; d is the lower bound; 
#	#'	scal is the slope of the linear middle part of the line; 
#	#'	xmid is x value at which y equals to (a-d)/2; g is the asymmetric factor.\cr
#	#'	Another assumption is that the 5pl lines for different
#	#'	batches should follow a similar pattern, but shifted horizontally left or
#	#'	right due to factors such as the changes in the standard reagents. In terms of the
#	#'	5pl function parameters, the 5pl lines for different batches are identical
#	#'	in a, d, scal and g, but different in xmid; That is to say that if we shift these 
#	#'	lines horizontally with the correct amounts, the lines can merge into one line.
#	#'	Therefore, we fit these lines together to one 5paremeter logistic function with
#	#'	identical a, d, scal and g, as well as one xmid for the reference line and one
#	#'	shift for each other lines. In this function, we use a linear model to generate
#	#'	the initial values for these shifts parameter of the fitting. \cr
#	#'	The algorithm takes a local linear model to generate the initial values like this,
#	#'	
#	#'	\itemize{
#	#'	\item {1. carefully choose the data set as the reference\cr}
#	#'		{
#	#'		all the lines shifted towards this one.
#	#'		the reference line could be anywhere in the data.
#	#'		Theoretically, it can
#	#'		be any one, but empirically we pick the one line, whose maximum value is closest
#	#'		to d/2 (65535/2). This way the reference line has the best coverage over the range of 
#	#'		5pl and easy to converge for the fitting.
#	#'		}
#	#'	\item {2. determine the initial k values\cr}
#	#'		{
#	#'		we simply compare the biggest Y in 
#	#'		in each data series. Two cases are possible \cr
#	#'			\itemize{
#	#'				\item {1}{ Ymax_i < Ymax_r, nothing to do }
#	#'				\item {2}{ Ymax_i >Ymax_r, find the max Yj_i in the series to be smaller then Ymax_r}
#	#'			}
#	#'		Now we have a pair (Yj_i, Xj_i).
#	#'		}
#	#'	\item {3.Determine where does this pair belongs to inside the Yr data line\cr}
#	#'		{
#	#'		Find (Y(k-1)_r, Yk_r), where Y(k-1)_r<Yj_i<Yk_r. Then 
#	#'		simply using a linear model to determine Xj_i_pred for Yj_i according
#	#'		to reference data . It could possible happen that Yj_i is not inside
#	#'		reference range, meaning Yj_i<Ymin_r. In this case, we simply assuming
#	#'		Yj_i == Ymin_r, using X(Ymin_r) as the predication and shift the data.
#	#'		Hopefully, there will not be many cases like this.
#	#'		}
#	#'}
#	#'To estimate the shift k, we take 
#	#\eqn{k=log(Xj_i/Xj_i_pre)}\cr
#	#'   k=log(Xj_i/Xj_i_pre)    \cr
#	#'Another note is that we specify as an input which batch to use as the "reference".
#	#'Then within the batch, we pick as indicated above the reference standard line
#	#'. The return value records the batch, the run and the plate as output. Also
#	#'the inits output has a number of elements identical to the total number of
#	#'plates.
#	#The accessary function for preparint intial values for nlsLM
#	#this is necessary, because we might have many data series.
#	#we need this for actually prepare the initial values
#	
#	#take in the input, before the prepareInput 
#	#and make the initial values in order to do nlsLM fitting
#	#it returns an array which has length identifical to the total series (row lengths)
#	#
	#'@param batches list of elisa_batch data 
	#'		
	#'@param ref.batch numeric the index of the reference batch. It is 1
	#'		by default.
	#'
	#'@return a data list contain the following elements,
	#' \itemize{
	#'	\item {inits, the initial values for the standard curves of all the plates\cr}
	#'		
	#'	\item {ref.ibatch, the index of the reference batch\cr}
	#'		{
	#'			This one is specified by the input ref.batch.
	#'		}
	#'	\item {ref.irun,the index of the reference run\cr}
	#'
	#'	\item {ref.index,the index of the reference line in the order
	#'			of the inits vector\cr}
	#'
	#'}
	# #'@export
	prepareInitsLM<-function(batches, ref.batch=1)
	{
		if(missing(batches))
		{
			stop("ERROR:please specify the input")
		}
		
		rbatch<-batches[[ref.batch]]
		#get reference batch need to find the reference line
		middles<-findMiddle.batch(rbatch)
		ref.middle<-middles$middle #middle point 
		
		inits<-c();
		#now go through every single line to get inits
		shift.indx<-0
		count<-0
		for(i in 1:length(batches))
		{
			ebatch<-batches[[i]]
			range.batch<-ebatch@range.ODs 
			OD.middle<-(range.batch[1]+range.batch[2])/2
			for(j in 1: ebatch@num.runs)
			{
				erun<-ebatch@runs[[j]]
				
				for(k in 1:erun@num.plates)
				{
					count<-count+1
					eplate<-erun@plates[[k]];
					xmiddle<- findMiddle.plate(eplate, OD.middle)
					if(middles$ind.run==j&&middles$ind.plate==k&&i==ref.batch)
					{
						inits<-c(inits,0)
						shift.index<-count;
					} else {
						inits<-c(inits,log(ref.middle/xmiddle))
					}						
				}
			}
		}
		
		list("inits"=inits, "ref.iRun"=middles$ind.run, 
			"ref.ibatch"=ref.batch, "ref.iplate"=middles$ind.plate,ref.index=shift.index)
	}

#find the middle point x/conc value according to the OD.middle
#using the simple linear model.	
findMiddle.plate<-function(eplate, OD.middle)
{
	m<-findUpLow.middle(eplate@data.std$OD, eplate@data.std$conc,OD.middle)
	middle<- reverseLookupX.LM(m$low,m$high,OD.middle)
}


#'@title Get the OD ranges (min/max)
#'@description Going through the list of  batches to get the OD range (min and max) 
#'@param batches list of batches data
#'@export
rangeOD<-function(batches)
{
	if(missing(batches))
	{
		stop("please specify the input ")
	}
	if(class(batches)!="list")
	{
		stop("please specify the correct input")
	}
	min.OD<-1000
	max.OD<- -1
	for(i in 1:length(batches))
	{
		range.ODs<-batches[[i]]@range.ODs
		if(min.OD< range.ODs[1])
		{
			min.OD<-range.ODs[1]
		}
		if(max.OD>range.ODs[2])
		{
			max.OD<-range.ODs[2]
		}
	}
	return (c(min.OD, max.OD))
}#

#summary and get one representative curve line 
#for each batch. But how??
#1)pool the x/concs and then get unique concs
#2)for conc, we obtain the average y/OD values
#we can do this by fitting one good 5pl line,
#but this might not be necessary. 
#so do it in a crude way, simply find one representative/middle
#line in each batch. 
#so we don't have to do it very accurately.
#in each batch
# #'@export
findMiddle.batch<-function(ebatch)
{
	rlen<-length(ebatch@runs)
	OD.min<-ebatch@range.ODs[1]
	OD.max<-ebatch@range.ODs[2]
	OD.middle<-(OD.max-OD.min)/2+OD.min;
	#now to go through each individual line
	#reverse calculate x coordinate of the middle point of 
	#each line xmid and then get the middle 
	#line as the representative of the batch
	#here we don't do 5pl fit, but just
	# a linear fit between the two points around
	# the middle point
	
	middles<-c();
	ind.run<-c();
	ind.plate<-c(); #index run and plate together rember which plate is for which middle points
	ind.count<-0
	for(i in 1:rlen)
	{
		ind.count<-ind.count+1;
		erun<-ebatch@runs[[i]]
		plen<-length(erun@plates)
		for(j in 1:plen)
		{
			eplate<-erun@plates[[j]]
			#find the range
			m<-findUpLow.middle(eplate@data.std$OD, eplate@data.std$conc,OD.middle)
			middles<-c(middles, reverseLookupX.LM(m$low,m$high,OD.middle))
			ind.run<-c(ind.run, i)
			ind.plate<-c(ind.plate,j)
		}
	}#end of run,
	#now we have the middle points. we need to find the middle line 
	ix<-order(middles)[ceiling(length(middles)/2)]
	plate<-ebatch@runs[[ ind.run[ix] ]]@plates[[ ind.plate[ix] ]];
	
	return(list(OD=plate@data.std$OD,conc=plate@data.std$conc, middle=middles[ix] , 
			ind.run=ind.run[ix], ind.plate=ind.plate[ix]))
}#end of function findMiddel.batch

reverseLookupX.LM<-function(p1, p2,ymid)
{
	#do linear fitting
	#y<-ax+b
	a<-(p1[2]-p2[2])/(p1[1]-p2[1])
	b<-p1[2]-a*p1[1]
	
	return((ymid-b)/a)
}

#to find two adjacent points contain OD.middle point 
#input, standard line seris with y (OD) and x (conc)
#return two points low and high contains middle point 
#we will do mean y over x to look up.
findUpLow.middle<-function(y, x, OD.middle)
{
	#first get mean y over x,
	tempDF<-aggregate(y, by=list(group=x),FUN="mean")
	y<-tempDF$x;
	x<-tempDF$group;
	#now start doing the job
	ind<-order(x)
	x.sort<-x[ind] #increasing
	y.sort<-y[ind] #increasing
	
	#now look up the range
	low<-c()
	high<-c()
	found.it<-FALSE
	for(i in 1:length(x))
	{
		if(y.sort[i]>OD.middle) #starting from low to high, we find one higher than OD.middle point
		{
			found.it<-TRUE;
			if(i==1) #special case, even the first one is larger than OD.middle, so we have out of range points, so what? no problem
			{
				low<-c(x.sort[i],y.sort[i])
				high<-c(x.sort[i+1],y.sort[i+1])
			} else {#this is the correct/normal case, OD.middle is in the middle of x, y
				low<-c(x.sort[i-1],y.sort[i-1])
				high<-c(x.sort[i],y.sort[i])
			}
			break;
		}
	}
	if(!found.it)  #means all the y in the series is smaller than 
	{
		low<-c(x.sort[length(x)-1],y.sort[length(x)-1]);
		
		high<-c(x.sort[length(x)],y.sort[length(x)]);
	}
	
	return(list(low=low, high=high))
}

# #'save the fitted regression model into the data
# #'	here we save the model into batch data
# #'	we save the regression model into batch level
# #'	save k/shift at plate level.
# #'	we also summarize the batch correction factor in this function
# #'
# #'	input
# #'		regModel contains the shifts and fitted parameters. remember 
#      the shifts do not include the reference line. so we need to 
#		insert it.
# #'
# #'@export

saveRegressionModel<-function(batches, regModel, mode=c("fix.both","fix.low", "fix.high","fix.all"),
			ref.index=1, a, d, xmid, scal, g, model=c("5pl","4pl"))
{
	if(missing(batches))
	{
		stop("ERROR:please specify the input batches")
	}
	if(class(batches)!="list")
	{
		stop("ERROR:please specify the input as list")
	}
	if(class(regModel)!="nls.lm")
	{
		stop("ERROR:please specify the input of the regression model")
	}
	model<-match.arg(model);
	ind<-0;
	pars<-c();
	ks<-c();
	mode<-match.arg(mode);
	if(model=="5pl")
	{
		switch(mode,
		"fix.both"={ind<-4; pars<-c(a, d,regModel$par[1:3]);ks<-regModel$par[-c(1:3)]},
		"fix.low"={ind<-5;pars<-c(a, regModel$par[1:4]);ks<-regModel$par[-c(1:4)]},
		"fix.high"={ind<-3;pars<-c(regModel$par[1], d,regModel$par[2:4]);ks<-regModel$par[-c(1:4)]},
		"fix.all"={ind<-1;pars<-c(a, d,xmid, scal, g);ks<-regModel$par}
		)
	} else { "4pl"
		g<-1;
		switch(mode,
		"fix.both"={ind<-3; pars<-c(a, d,regModel$par[1:2],g);ks<-regModel$par[-c(1:2)]},
		"fix.low"={ind<-4;pars<-c(a, regModel$par[1:3],g);ks<-regModel$par[-c(1:3)]},
		"fix.high"={ind<-2;pars<-c(regModel$par[1], d,regModel$par[2:3],g);ks<-regModel$par[-c(1:3)]},
		"fix.all"={ind<-1;pars<-c(a, d,xmid, scal, g);ks<-regModel$par}
		)
	}
	names(pars)<-c("a","d","xmid","scal","g")
	#now, let's insert zero into the array 
	if(ref.index>length(ks)+1||ref.index<=0)
	{
		stop("ERROR:shift.index out of range")
	}
	if(ref.index==1){
		ks<-c(0,ks)
	} else {
		if(ref.index==length(ks)+1){
			ks<-c(ks,0)
		} else {
			#insert in the middle
			ks<-c(ks[1:(ref.index-1)],0, ks[-(1:(ref.index-1))])
		}
	}
	#now we have everything ready, save data into plates
	count<-0;
	for(i in 1:length(batches))
	{
		normFactors<-c();
		for(j in 1:batches[[i]]@num.runs)
		{
			for(k in 1:batches[[i]]@runs[[j]]@num.plates)
			{
				count<-count+1;
				normFactors<-c(normFactors,ks[count])
				batches[[i]]@runs[[j]]@plates[[k]]@normFactor<-ks[count]
			}
		}
		batches[[i]]@pars<-pars;
		batches[[i]]@model.name<-model;
		#batches[[i]]@model.fit<-regModel;
		batches[[i]]@normFactor<-mean(normFactors);
	}
	
	return(batches)
}	

#align the standard by shifting toward reference and then plot
#then together to QC the fitting.
#graph.file: is the file name for the graph to be plotted.
#'@title Plot all batch data together
#'@description Plot the batch data together for visualization. 
#'@details If the data has been analysed, a fitted line will be drawn too. If
#'	there are more than one batches, each batch will be plotted with different color
#'	and different synmbols. Different batches will also be shifted/adjusted based on their
#'	"S" factor, and one single fitted line (based on the "reference" batch) will be plotted.
#'
#'@param batches list of batch data objects either raw or analyzed data.
#'@param graph.file characters as the output graph file name. If specified, a 
#'	SVG (*.svg) graph will be saved to the disk. Otherwise, the graph 
#'	will be send to the stdout.
#'
#'@return characters which specify the graph file name, if graph.file is specified. NULL
#'	otherwise.
#'
#'@examples
#' #load the library
#' library(ELISAtools)
#' 
#' #get file folder
#' dir_file<-system.file("extdata", package="ELISAtools")
#'
#' #load the data
#' batches<-loadData(file.path(dir_file,"design.txt"))
#'
#' #plot the raw batch data together
#' plotAlignData(batches);
#'
#'
#'@export		
plotAlignData<-function(batches, graph.file=NULL)
{
	if(missing(batches))
	{
		stop("please specify the input")
	}
	if(class(batches)!="list")
	{
		stop("the input should be a list")
	}
	if(!is.null(graph.file)&&graph.file!="")
	{
		#svg(filename=graph.file)
        png(filename=graph.file)
	}
	#plot 
	if(length(batches)<=0)
	{
		#no data
		cat("no data in the list");
		plot(c(1,1),c(2,2), type="n", xlab="conc", ylab="OD");
		return(NULL);
	}
	
	#now, let's do ploting expectation first
	pars<-batches[[1]]@pars;
	log_conc<-log(avoidZero(batches[[1]]@runs[[1]]@plates[[1]]@data.std$conc))
	x_min<-min(log_conc)
	x_max<-max(log_conc)
	
	x<-seq(x_min, x_max*1.1,by=(x_max-x_min)/1000);
	ymin<-batches[[1]]@range.ODs[1]
	ymax<-batches[[1]]@range.ODs[2]
	flag.analyzed<-T;
	if((batches[[1]]@model.name=="5pl"||batches[[1]]@model.name=="4pl")&&length(pars)==5)
	{
		#pars<-pars+c(0,0,-1*batches[[1]]@normFactor,0,0)
		y<-f5pl(pars, x)
		plot(x,y, type="l",xlab="conc", ylab="OD", lwd=2, 
			lty=2,col=1, main=paste0("ELISA Batch Data:","Fitted and Adjusted"));
		flag.analyzed<-T;
	} else { #for non-analyzed data.
		plot(c(x_min,x_max),c(ymin,ymax), type="n",xlab="conc", ylab="OD", lwd=2, 
			lty=2,col=1, main=paste0("ELISA Batch Data:","Raw"));
		flag.analyzed<-F;
	}
	#y<-f5pl(pars, x)
	#plot(x,y, type="l",xlab="conc", ylab="OD", lwd=2, 
	#		lty=2,col=1, main="ELISA Aligned Data");
	
	#plot lines
	#count<-0;
	batchIDs<-c();
	for(i in 1:length(batches))
	{
		for(j in 1:batches[[i]]@num.runs)
		{
			for(k in 1:batches[[i]]@runs[[j]]@num.plates)
			{
				#count<-count+1;
				conc<-log(avoidZero(batches[[i]]@runs[[j]]@plates[[k]]@data.std$conc))
				fac<-0;
				if(flag.analyzed&&!is.na(batches[[i]]@runs[[j]]@plates[[k]]@normFactor))
				{
					fac<-batches[[i]]@runs[[j]]@plates[[k]]@normFactor
				}
				conc<-conc+fac;
				points(conc, batches[[i]]@runs[[j]]@plates[[k]]@data.std$OD, col=i+1, pch=(i+1)%%25)
				dF<-aggregate(batches[[i]]@runs[[j]]@plates[[k]]@data.std$OD,
						by=list(conc=batches[[i]]@runs[[j]]@plates[[k]]@data.std$conc),FUN=mean)
				lines(log(avoidZero(dF$conc))+fac,dF$x,col=i+1,lty=3);
			}
		}
		batchIDs<-c(batchIDs,batches[[i]]@batchID)
	}
	batchIDs<-c("fitted line",batchIDs)
	legend(x_min,ymax/2,#pars[2]/2,
			batchIDs, lty=c(2, rep(3,length(batchIDs)-1)),col=c(1:length(batchIDs)), 
			pch=c(-1,(1:length(batchIDs))+1), lwd=c(2,rep(1,length(batchIDs)-1))
			)
	if(!is.null(graph.file)&&graph.file!="")
	{
		dev.off();
	}
	return(graph.file);
}

#do not align the standards. simply plot the data by batch. 
#no shifting
#graph.file is the file name for the graph to be plotted.
#'@title Plot ELISA data for one batch
#'@description Plot the individual batch data for visualization. 
#'@details If the data has been analysed, a fitted line will be drawn too. 
#'
#'@param batch batch data objects with either raw or analyzed data.
#'@param graph.file characters as the output graph file name. If specified, a 
#'	SVG (*.svg) graph will be saved to the disk. Otherwise, the graph 
#'	will be send to the stdout.
#'
#'@return characters which is the graph file name, if graph.file is specified. NULL
#'	otherwise.
#'
#'@examples
#' #load the library
#' library(ELISAtools)
#' 
#' #get file folder
#' dir_file<-system.file("extdata", package="ELISAtools")
#'
#' #load the data
#' batches<-loadData(file.path(dir_file,"design.txt"))
#'
#' #plot the raw batch 1 data
#' plotBatchData(batches[[1]]);
#'
#'
#'@export
plotBatchData<-function(batch, graph.file=NULL)
{
	if(missing(batch))
	{
		stop("ERROR: please specify the batch input")
	}
	
	if(class(batch)!="elisa_batch")
	{
		stop("ERROR: please specify the input as elisa_batch data")
	}
	
	pars<-batch@pars;
	log_conc<-log(avoidZero(batch@runs[[1]]@plates[[1]]@data.std$conc))
	x_min<-min(log_conc)
	x_max<-max(log_conc)
	
	x<-seq(x_min, x_max*1.1,by=(x_max-x_min)/1000);
	ymin<-batch@range.ODs[1]
	ymax<-batch@range.ODs[2]
	
	if(!is.null(graph.file)&&graph.file!="")
	{
		#svg(filename=graph.file)
        png(filename=graph.file)
	}
	
	if((batch@model.name=="5pl"||batch@model.name=="4pl")&&length(pars)==5)
	{
		pars<-pars+c(0,0,-1*batch@normFactor,0,0)
		y<-f5pl(pars, x)
		plot(x,y, type="l",xlab="conc", ylab="OD", lwd=2, 
			lty=2,col=1, main=paste0("ELISA Batch Data:",batch@batchID));
	} else { #for non-analyzed data.
		plot(c(x_min,x_max),c(ymin,ymax), type="n",xlab="conc", ylab="OD", lwd=2, 
			lty=2,col=1, main=paste0("ELISA Batch Data:",batch@batchID));
	}
	
	count<-0;
	runID<-c();
	for(j in 1:batch@num.runs)
	{
		for(k in 1:batch@runs[[j]]@num.plates)
		{
			count<-count+1;
			conc<-log(avoidZero(batch@runs[[j]]@plates[[k]]@data.std$conc))
			#fac<-batch@runs[[j]]@plates[[k]]@normFactor
			#conc<-conc+fac;
			points(conc, batch@runs[[j]]@plates[[k]]@data.std$OD, col=j+1, pch=count%%25)
			dF<-aggregate(batch@runs[[j]]@plates[[k]]@data.std$OD,
					by=list(conc=batch@runs[[j]]@plates[[k]]@data.std$conc),FUN=mean)
			lines(log(avoidZero(dF$conc)),dF$x,col=j+1,lty=3);
		}
		runID<-c(runID, j)
	}
	runID<-c("fitted batch mean", paste0("Run_",runID));
	#batchIDs<-c(batchIDs,batches[[i]]@batchID)
	legend(x_min,ymax/2,#pars[2]/2,
			runID, lty=c(2, rep(3,length(runID)-1)),col=c(1:length(runID)), 
			pch=rep(-1,length(runID)), lwd=c(2,rep(1,length(runID)-1))
			)
	if(!is.null(graph.file)&&graph.file!="")
	{
		dev.off();
	}
	return(graph.file);
}
###
#write html report
#'@title Report ELISA data in HTML format.
#'@description Writting the ELISA analysis results by batch in HTML format. 
#'@param batches list of elisa batch data objects. The data can be raw or after 
#'	analyzed and batch-corrected.
#'@param file.name character string denoting the report file. The file will be
#'		written in HTML format.
#'@param file.dir character string denoting the directory to save the report. 
#'@param desc character string describing the project and experiment. Will be 
#'		written into the report.
#'@return the function returns NULL. But it will save the html report to the disk.
#' Therefore, it is IMPORTANT to specify a directory you have write permission to
#'	run this function. 
#'
#'@examples
#'#R code to run 5-parameter logistic regression on ELISA data
#'#load the library
#'library(ELISAtools)
#'
#'##
#'#get file folder
#'dir_file<-system.file("extdata", package="ELISAtools")
#'
#'batches<-loadData(file.path(dir_file,"design.txt"))
#'
#'
#'#----IMPORTANT-----
#'#please make sure you have the write permission to save the html report
#'reportHtml(batches,file.dir=tempdir());
#'
#'@seealso  \code{\link{elisa_batch}} \code{\link{elisa_run}}
#'			\code{\link{elisa_plate}} 
#'
#'@export 
#assuming all the analyses like reading, fitting and predicting have been accomplished 
#updated 8/24/2018, now it could take care of raw data without analysis
reportHtml<-function(batches, file.name="report", file.dir=".", desc="")
{
	if(missing(batches))
	{
		stop("ERROR:please specify the input data")
	}
	
	if(class(batches)!="list")
	{
		stop("ERROR:please specify the input data as a list of elisa_batch objects")
	}
	
	if(!file.exists(file.dir))
	{
		stop("ERROR:the specified directory doesn't exist. please check")
	}
		
	if(file.exists(file.path(file.dir, paste0(file.name,"html"))))
	{
		cat("The specified file exists and will be overwritten");
	}
	file.dir<-file.path(file.dir);
		#now start 
	HTMLStart(outdir=file.dir, filename=file.name,
		extension="html", echo=F, HTMLframe=FALSE#, append=TRUE
		)
	x<-HTML.title("ELISA Data Analysis Tool Report", HR=1);
	x<-{
		if(nchar(desc)==0){
			desc<-paste0("ELISA batch analysis done on ",date());
		}
	};
	x<-HTML.title(desc, HR=3);
	#summary(pars)
	#writing fitting information and QC first
	x<-{
		if(length(batches)==0)
		{
			data.frame("empty input data found and no analysis reported")
			HTMLStop();
			return(NULL);
		}
	};
	x<-HTML(data.frame(desc="Regression Model:", content=batches[[1]]@model.name));
	if(!is.null(batches[[1]]@model.name)&&batches[[1]]@pars!=-1){
		x<-HTML("Model Parameters:");
		y<-{
			if(batches[[1]]@model.name=="5pl"){
				#x<-HTML("Equation 1:a+(d-a)/(1+exp((xmid-x)/scal))^g")
				x<-HTML(batches[[1]]@pars,caption="Equatin 1:y=a+(d-a)/(1+exp((xmid-x)/scal))^g",
						captionalign="top",nsmall=2)
			} else { #"4pl" in this case
				#x<-HTML()
				x<-HTML(batches[[1]]@pars[1:4],caption="Equation 1:y=a+(d-a)/(1+exp((xmid-x)/scal))",
						captionalign="top",nsmall=2)
			}
		};
		###now adding the parameters for regular scale
		y<-{
			if(batches[[1]]@model.name=="5pl"){
				param_2ndform<-c(batches[[1]]@pars[c("a","d")],exp(batches[[1]]@pars["xmid"]),-1/batches[[1]]@pars["scal"],batches[[1]]@pars["g"]);
				names(param_2ndform)<-c("D","A","C","B","g");
				#x<-HTML()
				x<-HTML(param_2ndform,caption="Equation 1:y=D+(A-D)/(1+(x/C)^B)^g",
						captionalign="top",nsmall=2);
			} else { #"4pl" in this case
				param_2ndform<-c(batches[[1]]@pars[c("a","d")],exp(batches[[1]]@pars["xmid"]),-1/batches[[1]]@pars["scal"]);
				names(param_2ndform)[1:4]<-c("D","A","C","B");
				#x<-HTML(, )
				x<-HTML(param_2ndform,caption="Equation 2: y=D+(A-D)/(1+(x/C)^B)",
						captionalign="top",nsmall=2)
				
			}
		};
		x<-HTML("S Factors:");
		#make one data frame to output the s factors
		bIds<-c();
		bNfac<-c();
		for(j in 1:length(batches)){
			batch<-batches[[j]];
			bIds<-c(bIds,batch@batchID);
			bNfac<-c(bNfac,batch@normFactor);
		}
			#dfm<-data.frame();
		x<-HTML(data.frame("Batch"=bIds,"S_Factor"=bNfac));
	}
	x<-HTMLhr();
	x<-HTML.title("Model fitting QC", HR=3);
	#if(
	fname.suffix<-as.numeric(format(Sys.time(), "%OS3"))*1000 
	#graphName<-file.path(file.dir,paste0("aligned_",fname.suffix,".svg"));
    graphName<-file.path(file.dir,paste0("aligned_",fname.suffix,".png"));
	x<-plotAlignData(batches, graphName);
	x<-HTMLInsertGraph(graphName,file=file.path(file.dir,paste0(file.name,".html")));
	#x<-HTMLplot() ;
	x<-HTMLhr();
	x<-Sys.sleep(0.95);
	x<-HTML.title("Data Output:", HR=3);
	for(i in 1:length(batches))
	{
		batch<-batches[[i]];
		
		x<-HTML.title(paste0("Batch:",batch@batchID,"; S Factor:", format(batch@normFactor,digit=3,nsmall=2)), HR=2);
		#graphName<-file.path(file.dir,paste0("batch_",i,"_",fname.suffix, ".svg"));
        graphName<-file.path(file.dir,paste0("batch_",i,"_",fname.suffix, ".png"));
		x<-plotBatchData(batch, graphName);
		#cat("graph name is:", graphName,"\n");
		x<-HTMLInsertGraph(graphName,file=file.path(file.dir,paste0(file.name,".html")));
		#cat("file name is :", file.path(file.dir,paste0(file.name,".html")),"\n");
		#x<-HTMLplot() ;
		#write data.
		x<-Sys.sleep(0.1);
		#x<-HTMLhr();
		for(j in 1:batch@num.runs)
		{
			x<-HTML.title(paste0("Run #",j,"\t",batch@runs[[j]]@date), HR=3);
			x<-HTMLhr();
			for(k in 1:batch@runs[[j]]@num.plates)
			{
				x<-HTML(paste0("RUN_# ",j,"\t",batch@runs[[j]]@date, "\tplate_#",k,";\tS Factor:",format(batch@runs[[j]]@plates[[k]]@normFactor,digit=3,nsmall=2)));
				if(!is.null(batch@runs[[j]]@plates[[k]]@data.unknown)&&dim(batch@runs[[j]]@plates[[k]]@data.unknown)[1]!=0)
				{
					x<-HTML(batch@runs[[j]]@plates[[k]]@data.unknown, nsmall=3, 
							caption="Sample of Unknown concentration", captionalign="top"
							)
				}
				if(!is.null(batch@runs[[j]]@plates[[k]]@mdata.unknown)&&dim(batch@runs[[j]]@plates[[k]]@mdata.unknown)[1]!=0)
				{
					x<-HTML(batch@runs[[j]]@plates[[k]]@mdata.unknown, nsmall=3, 
						caption="Sample of Unknown concentration(Mean)", captionalign="top")
				}
				x<-HTMLhr();
			}
		}
		
	}
	HTMLStop();
	fn<-file.path(file.dir, paste0(file.name,".txt"));
	saveDataText(batches, fn);
	cat("\nAn html reprot,\"",paste0(file.name,".html\", has been generate.\n"));
	cat("\nA text file,\"", paste0(file.name, ".txt\", has been gerated.\n")); 
	cat("\n");
	#return(NULL);
}