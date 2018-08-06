
#'@import R2HTML

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
#'@title prepare the input for regressoin
#'@description  blank for now
#'@details  blank for now
#'@param batches list contain the ELISA data, arranged in 
#'	batches. Each element of list contains a batch (list)
#'	data, and each batch contains multiple elisa_run 
#'	objects \code{\link{elisa_plate}}
#'	
#'@param ref.ID character the reference batch ID. all other 
#'	batches are shifted/corrected towards this reference batch
#' 
#'@return list of parameters that will feed in to do regression.
#'@seealso \code{\link{elisa_plate}}
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

#'@title prepare initial values for fitting shifts 
	#'@description generate the initial values for fitting shifts with 
	#'	a model of the 5-parameter logistic function.
	#'@details this is a more complicated way to prepare the initials for shifting. 
	#'	The assumption is that the signal strengths of a feature in protein microarray
	#'	(ProtoArray, invitrogen) acquired with different PMT gains follow 
	#'	a 5-parameter logistic function (5pl) as \cr
	#'	\ifelse{html}{\out{<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML"></script>
	#'						$$y = a + {{d-a} \over {(1 + e^ {{(xmid - x)} \over scal})}^g} $$
	#'				 	}}{\deqn{y=a+\frac{d-a}{(1+e^(\frac{xmid-x}{scal}))^g}}{ASCII}}
	#'	\cr 
	#'	where a is the upper bound; d is the lower bound); 
	#'	scal is the slope of the linear middle part of the line; 
	#'	xmid is x value at which y equals to (a-d)/2; g is the asymmetric factor.\cr
	#'	On the same array, the 5pl lines for different
	#'	feature proteins should follow a similar pattern, but shifted horizontally left or
	#'	right based on the level of the secondary antibodies bound, which is determined
	#'	by the level of printed feature proteins as well as the interaction between
	#'	the feature protein and the testing antibody/protein. In terms of the
	#'	5pl function parameters, the 5pl lines for different feature proteins are identical
	#'	in a, d, scal and g, but different in xmid; That is to say that if we shift these 
	#'	lines horizontally with the correct amounts, the lines can merge into one line.
	#'	Therefore, we fit these lines together to one 5paremeter logistic function with
	#'	identical a, d, scal and g, as well as one xmid for one reference line and one
	#'	shift for each other lines. In this function, we use a linear model to generate
	#'	the initial values for these shifts parameter of the fitting. \cr
	#'	The algorithm takes a local linear model to generate the initial values like this,
	#'	
	#'	\itemize{
	#'	\item {1. carefully choose the data set as the reference\cr}
	#'		{
	#'		all the lines shifted towards this one.
	#'		the reference line could be anywhere in the data.
	#'		Theoretically, it can
	#'		be any one, but empirically we pick the one line, whose maximum value is closest
	#'		to d/2 (65535/2). This way the reference line has the best coverage over the range of 
	#'		5pl and easy to converge for the fitting.
	#'		}
	#'	\item {2. determine the initial k values\cr}
	#'		{
	#'		we simply compare the biggest Y in 
	#'		in each data series. Two cases are possible \cr
	#'			\itemize{
	#'				\item {1}{ Ymax_i < Ymax_r, nothing to do }
	#'				\item {2}{ Ymax_i >Ymax_r, find the max Yj_i in the series to be smaller then Ymax_r}
	#'			}
	#'		Now we have a pair (Yj_i, Xj_i).
	#'		}
	#'	\item {3.Determine where does this pair belongs to inside the Yr data line\cr}
	#'		{
	#'		Find (Y(k-1)_r, Yk_r), where Y(k-1)_r<Yj_i<Yk_r. Then 
	#'		simply using a linear model to determine Xj_i_pred for Yj_i according
	#'		to reference data . It could possible happen that Yj_i is not inside
	#'		reference range, meaning Yj_i<Ymin_r. In this case, we simply assuming
	#'		Yj_i == Ymin_r, using X(Ymin_r) as the predication and shift the data.
	#'		Hopefully, there will not be many cases like this.
	#'		}
	#'}
	#'to estimate the shift k, we take 
	#'\eqn{k=log(Xj_i/Xj_i_pre)}
	
	#The accessary function for preparint intial values for nlsLM
	#this is necessary, because we might have many data series.
	#we need this for actually prepare the initial values
	
	#take in the input, before the prepareInput 
	#and make the initial values in order to do nlsLM fitting
	#it returns an array which has length identifical to the total series (row lengths)
	#
	#'@param y dataframe this is raw data and has NOT been prepared. 
	#'		it has row as genes and col as pmts
	#'@param x numeric vector It contains the PMT settings for each array. 
	#'@param data.aggregated Parameter indicates whether the data has been 
	#'		aggregated, mean the two repeats has been averaged. 
	#'		By default it is not. In this case, we will generate 
	#'		the two repeat the identical shift, k.
	#
	#'@param ref.line list mean to use an external line as the reference line.
	#'		if missing then, use the internal line as the reference.
	#'		It is a list containing two items. The first one is the ref.index 
	#'		for the index of the reference line in the first data set. The
	#'		second is the actually lines. this could be two line for unaggregated
	#'		data or one for aggregated data.
	#'
	#'@return a data list contain two things. First is the ref.index, the index
	#'		of the line in the data set. If we are using the external reference 
	#'		line, then this is -1. The second is estimated initial values for
	#'		shifts. The shifts INCLUDE the shift for the reference line with
	#'		a value of zero. 
	#'		Note: 1) the prepared inits vector always is in a format of aggregated data
	#'		; 2) in log'ed fold; and
	#'		the output doesn't include ref.line when we reference to the external lines
	#'@export
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
##'@export
#get X Y data from elisa plate data
##'@title to obtain X (conc) and Y (OD) from the elisa_plate data
##'@description obtain standard X and Y data to put them into 
##'	vector, so that to make them  
#getXY.std<-function(plate)
#{
#		return list(x=plate@data.std$OD, y=plate@data.std$conc)
#}

#getXYs.std<-function(eplates)
#{
	
#	for(
#}
#go through the runs in the batch to get 
#the range of ODs to make a guess of a and d
#of the 5pl 
#'@export
rangeOD<-function(batches)
{
	if(missing(batch))
	{
		stop("please specify the input ")
	}
	if(class(ebatch)!="list")
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
#'@export
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
	plate<-batch@runs[[ ind.run[ix] ]]@plates[[ ind.plate[ix] ]];
	
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

##'save the fitted regression model into the data
##'	here we save the model into batch data
##'	we save the regression model into batch level
##'	save k/shift at plate level.
##'	we also summarize the batch correction factor in this function
##'
##'	input
##'		regModel contains the shifts and fitted parameters. remember 
#      the shifts do not include the reference line. so we need to 
#		insert it.
#'@export

saveRegressionModel<-function(batches, regModel, mode=c("fix.both","fix.low", "fix.high","fix.all"),
			ref.index=1, a, d, xmid, scal, g)
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
	ind<-0;
	pars<-c();
	ks<-c();
	mode<-match.arg(mode);
	switch(mode,
	"fix.both"={ind<-4; pars<-c(a, d,regModel$par[1:3]);ks<-regModel$par[-c(1:3)]},
	"fix.low"={ind<-5;pars<-c(a, regModel$par[1:4]);ks<-regModel$par[-c(1:4)]},
	"fix.high"={ind<-3;pars<-c(regModel$par[1], d,regModel$par[2:4]);ks<-regModel$par[-c(1:4)]},
	"fix.all"={ind<-1;pars<-c(a, d,xmid, scal, g);ks<-regModel$par}
	)
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
		batches[[i]]@model.name<-"5pl";
		#batches[[i]]@model.fit<-regModel;
		batches[[i]]@normFactor<-mean(normFactors);
	}
	
	return(batches)
}	

#align the standard by shifting toward reference and then plot
#then together to QC the fitting.
#'@export		
plotAlignData<-function(batches)
{
	if(missing(batches))
	{
		stop("please specify the input")
	}
	if(class(batches)!="list")
	{
		stop("the input should be a list")
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
	y<-f5pl(pars, x)
	plot(x,y, type="l",xlab="conc", ylab="OD", lwd=2, 
			lty=2,col=1, main="ELISA Aligned Data");
	
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
				fac<-batches[[i]]@runs[[j]]@plates[[k]]@normFactor
				conc<-conc+fac;
				points(conc, batches[[i]]@runs[[j]]@plates[[k]]@data.std$OD, col=i+1, pch=i+1)
				dF<-aggregate(batches[[i]]@runs[[j]]@plates[[k]]@data.std$OD,
						by=list(conc=batches[[i]]@runs[[j]]@plates[[k]]@data.std$conc),FUN=mean)
				lines(log(avoidZero(dF$conc))+fac,dF$x,col=i+1,lty=3);
			}
		}
		batchIDs<-c(batchIDs,batches[[i]]@batchID)
	}
	batchIDs<-c("fitted line",batchIDs)
	legend(x_min,pars[2]/2,batchIDs, lty=c(2, rep(3,length(batchIDs)-1)),col=c(1:length(batchIDs)), 
			pch=c(-1,(1:length(batchIDs))+1), lwd=c(2,rep(1,length(batchIDs)-1))
			)
}

#do not align the standards. simply plot the data by batch. 
#no shifting
#'@export
plotBatchData<-function(batch)
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
	pars<-pars+c(0,0,-1*batch@normFactor,0,0)
	y<-f5pl(pars, x)
	plot(x,y, type="l",xlab="conc", ylab="OD", lwd=2, 
			lty=2,col=1, main=paste0("ELISA Batch Data:",batch@batchID));
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
			points(conc, batch@runs[[j]]@plates[[k]]@data.std$OD, col=j+1, pch=count)
			dF<-aggregate(batch@runs[[j]]@plates[[k]]@data.std$OD,
					by=list(conc=batch@runs[[j]]@plates[[k]]@data.std$conc),FUN=mean)
			lines(log(avoidZero(dF$conc)),dF$x,col=j+1,lty=3);
		}
		runID<-c(runID, j)
	}
	runID<-c("fitted batch mean", paste0("Run_",runID));
	#batchIDs<-c(batchIDs,batches[[i]]@batchID)
	legend(x_min,pars[2]/2,runID, lty=c(2, rep(3,length(runID)-1)),col=c(1:length(runID)), 
			pch=rep(-1,length(runID)), lwd=c(2,rep(1,length(runID)-1))
			)
}
###
#write html report

#'@export
#assuming all the analyses like reading, fitting and predicting have been accomplished 
reportHtml<-function(batches, file.name="report.html", file.dir=".", desc="")
{
	if(missing(batches))
	{
		stop("ERROR:please specify the input data")
	}
	
	if(class(batches)!="list")
	{
		stop("ERROR:please specify the input data as a list of elisa_batch objects")
	}
	
	if(!file.exists("."))
	{
		stop("ERROR:the specified directory doesn't exist. please check")
	}
	
	
	if(file.exists(file.path(file.dir, file.name)))
	{
		cat("The specified file exists and will be overwritten");
	}
	
		#now start 
	HTMLStart(outdir=file.dir, file=file.name,
		extension="html", echo=F, HTMLframe=FALSE#, append=TRUE
		)
	x<-HTML.title("ELISA Data Analysis Tool Report", HR=1);
	x<-{if(nchar(desc)==0){
		desc<-paste0("ELISA batch analysis done on ",date());
	}};
	x<-HTML.title(desc, HR=3);
	#summary(pars)
	#writing fitting information and QC first
	x<-{if(length(batches)==0)
	{
		data.frame("empty input data found and no analysis reported")
		HTMLStop();
		return(NULL);
	}};
	x<-HTML(data.frame(desc="Regression Model:", content=batches[[1]]@model.name));
	
	x<-HTML("Model Parameters:");
	x<-HTML(batches[[1]]@pars)
	x<-HTMLhr();
	x<-HTML.title("Model fitting QC", HR=3);
	x<-plotAlignData(batches);
	x<-HTMLplot() ;
	x<-HTMLhr();
	x<-Sys.sleep(0.6);
	x<-HTML.title("Data Output:", HR=3);
	for(i in 1:length(batches))
	{
		batch<-batches[[i]];
		
		x<-HTML.title(paste0(batch@batchID,":"), HR=2);
		x<-plotBatchData(batch);
		x<-HTMLplot() ;
		#write data.
		x<-Sys.sleep(0.9);
		#x<-HTMLhr();
		for(j in 1:batch@num.runs)
		{
			x<-HTML.title(paste0("Run_#",j), HR=3);
			x<-HTMLhr();
			for(k in 1:batch@runs[[j]]@num.plates)
			{
				x<-HTML(paste0("RUN_#",j, ":plate_#",k));
				x<-HTML(batch@runs[[j]]@plates[[k]]@data.unknown, nsmall=3, 
						caption="Sample of Unknown concentration", captionalign="top"
						)
				x<-HTML(batch@runs[[j]]@plates[[k]]@mdata.unknown, nsmall=3, 
						caption="Sample of Unknown concentration(Mean)", captionalign="top")
				x<-HTMLhr();
			}
		}
		
	}
	HTMLStop();
	
}