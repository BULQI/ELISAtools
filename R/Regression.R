#Regression module for doing ELISA analysis with batch correction.
#it is in this module to do the actually regression and prediction 
#	based on the regression model.
#In this model, we will calling on the functions in PAST for 
# fitting and shifting.


#import to use the fitting functions
# #'@import PAST
#'@import minpack.lm
# ##'@import stringi
#
#'@title Fit 5- or 4-parameter logistic function
#'@description Fit 5- or 4-parameter logistic function to estimate the parameters
#'	by pooling the standard curves from all batches
#'@details In this fitting, we first "guess" the initial values and then estimate
#'	the parameters based on 5- or 4-parameter function by shifting every single standard 
#'	curves towards the reference line. 
#'  We are reasoning that the intra-batch and 
#'	inter-batch factors affect the curve similarly by shifting the curve 
#'	left or right without changing its shapes. So we combine them together to
#'	to fit one single reference curve. To model the inter-batch effects, we
#'	take the average of the shifts of curves withine each batch, and use it
#'	to correct/normalize between different batches. \cr
#'	To summerize, each individual curve has its own shifts, which contains 
#'	the information about intra- and inter-batch effects. each batch has
#'	one batch level shifts (S Factor), which is an average of shifts of 
#'	curves within its batch and contains information about inter-batch effects.
#'	When we try to normalize between batches, we will apply the batch level shift to 
#'	all the curves within the same batch.
#'@param pars numeric vector initial values to estimate the paramters
#'@param a numeric the initial value for a (the lower limit of the function)  
#'@param d numeric the initial value for d (the upper limit of the function)
#'@param batches list of the batch data used to fit the model
#'@param refBatch.ID numeric or string indicating the 
#'	reference batch, by default is set to be the first one.
#'@param model characters to indicate either 5-parameter logistic
#'	function (5pl, default one) or 4-parameter logistic (4pl) to 
#'	be used in the fitting.
#'@return the batch data with the fitted model 
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
#'
#'#do fitting. model will be written into data set.
#'batches<-runFit(pars=pars,  batches=batches, refBatch.ID=1, model=model  )
#' 
#' @references Feng, et al 2018 \doi{10.1101/483800}
#'@export
#note 1) the outside caller need to prepare for data aggregation, in here
#we do data aggregation.
# 2) we do shift index here too.
#
runFit<-function(pars, a, d, batches, refBatch.ID=1,model=c("5pl","4pl")  )
{	#cat("\tinitial fitting for shifts and parameters.....\n")
	#inits, y, x,
	#prepare, y and x, 
	input<-prepareRegInput(batches);
	y<-input$y
	x<-input$x
	num.std<-input$num.std;
	#prepare the min and max as a and d, respectively
	if(missing(a))
	{
		a<-min(y)
	}
	if(missing(d))
	{
		d<-max(y)
	}
	model<-match.arg(model);
	if(model=="5pl")
	{
		if(length(pars)>3)
		{
			warnings("more parameters specified for 5-parameter logistic model and ignored\n");
		}
		if(length(pars)<3)
		{
			stop("not enought parameters specified for 5-parameter logisitic model, please check\n");
		}
	}
	if(model=="4pl")
	{
		if(length(pars)>2)
		{
			warnings("more parameters specified for 4-parameter logistic model and ignored\n");
		}
		if(length(pars)<2)
		{
			stop("not enought parameters specified for 4-parameter logistic model, please check\n");
		}
	}
	#now let's check for zero, since we need to do log
	x<-avoidZero(x)
	initsLM<-prepareInitsLM(batches, refBatch.ID)
	
	parStart<-c(pars,initsLM$inits[- initsLM$ref.index])
	fm<-NULL;
	if(model=="5pl"){
		fm<-nls.lm(par=parStart, 
				#fn=gainAdjust.fnRes4pShift,  #<----4pl, varying a
				fn=fnRes3pShift, #<----3pL, keeping a and d constant
				#y=log(yinput), 
				y=y,
				x=log(x),# xlen=xlen, aggregated=data.aggregated, ylog=data.ylog,
				a=a,d=d,
				shift.index=initsLM$ref.index,
				num.std=num.std,
				#model.weight="uniform", order=1, 
				control = nls.lm.control(nprint=1,maxiter=50)
			)
	} else { "4pl"
		fm<-nls.lm(par=parStart, 
				#fn=gainAdjust.fnRes4pShift,  #<----4pl, varying a
				fn=fnRes2pShift, #<----2paremeter for 4pL, keeping a and d and g constant
				#y=log(yinput), 
				y=y,
				x=log(x),# xlen=xlen, aggregated=data.aggregated, ylog=data.ylog,
				a=a,d=d,
				shift.index=initsLM$ref.index,
				num.std=num.std,
				#model.weight="uniform", order=1, 
				control = nls.lm.control(nprint=1,maxiter=50)
			)
	}	
	#now that we have regression model fitted, need to save the data into batches 
	#saveRegressionModel(batches, RegModel, mode, ref.index, a, b);
	batches<-saveRegressionModel(batches=batches,regModel=fm , a=a, d=d, ref.index=initsLM$ref.index,model=model)	
}

	##Residual function used by nlsLM it will estimate the
	## ks for shifting/aligning x input.
#pars=c(xmid, scal, g, k1, k2,...), notes: these input ks doesn't include the reference one, which always be zero.
#aggregate
# #'@export
fnRes3pShift<-function( pars, #parameters 
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		a, #par a, smallest y log
		d, #par d, largest y log
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		#g #the asymetric factor
		#xlen, #the length of distinct xs, x<-c(250,300,350,400,450,500,550,600,650), then xlens=9
		shift.index=1, #this is the index of the data set in the y array, to which all other data series shift.
		num.std #,this is num of data points in each standard, so we can use this information to expand the k (shifts) of each standard curve
		#ylog=TRUE,
		#aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		#, model.weight="sqrt", order=1 #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		
		
		#now also check for data consistency between inits and data
		#this is 3 parameter fitting,
		inits<-pars[-c(1:3)]
		#first insert zero into inits
		if(shift.index>length(inits)+1||shift.index<=0)
		{
			stop("ERROR:shift.index out of range")
		}
		if(shift.index==1){
			inits<-c(0,inits)
		} else {
			if(shift.index==length(inits)+1){
				inits<-c(inits,0)
			} else {
				#insert in the middle
				inits<-c(inits[1:(shift.index-1)],0, inits[-(1:(shift.index-1))])
			}
		}
		
		#now check the data integrity
		if(length(inits)!=length(num.std))
		{
			stop("ERROR: number of inits doesn't aggree with num.std, please check")
		}
		if(length(x)!=sum(num.std))
		{
			stop("ERROR: number of standards doesn't equal to num.std specified")
		}
		#now seems good the data input
		inits.expand<-rep(inits, num.std);
		xinput<-x+inits.expand
		
		xmid<-pars[1]
		scal<-pars[2]
		g<-pars[3]
		if(missing(a))
		{
			a<-min(y)
		}
		if(missing(d))
		{
			d<-max(y)
		}
		
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		#m<-weight.matrix(pred,model.weight,order)
		(y-pred)  #/m #sqrt(abs(pred))#^2#*(pred)
	}
	
	##here we use 5 parameter as 4pl by setting g==1;
	##Residual function used by nlsLM it will estimate the
	## ks for shifting/aligning x input.
#pars=c(xmid, scal, g, k1, k2,...), notes: these input ks doesn't include the reference one, which always be zero.
#aggregate
# #'@export
fnRes2pShift<-function( pars, #parameters 
		y, #the response, assuming log transformed
		x, #the independent variable, log-transformed
		#k, #vector of k to changing x 
		a, #par a, smallest y log
		d, #par d, largest y log
		#xmid, #x value at the middel point of y
		#scal, #slope at xmid
		#g, #the asymetric factor
		#xlen, #the length of distinct xs, x<-c(250,300,350,400,450,500,550,600,650), then xlens=9
		shift.index=1, #this is the index of the data set in the y array, to which all other data series shift.
		num.std #,this is num of data points in each standard, so we can use this information to expand the k (shifts) of each standard curve
		#ylog=TRUE,
		#aggregated=TRUE #indicated whether the data has been aggregated, mean the two repeats has been averaged
	#				by default it is not. So in this case we need to give the two repeat the identical shift, k
		#, model.weight="sqrt", order=1 #order is used for generation of exponential weight matrix, the order of the power. Not used if other type of weight matrix is applied.
		)
	{
		####check the input to make sure the input are in the correct format
		if(length(x)!=length(y))
		{
			stop("input x and y are not equal in length")
		}
		
		
		#now also check for data consistency between inits and data
		#this is 3 parameter fitting,
		inits<-pars[-c(1:2)]
		#first insert zero into inits
		if(shift.index>length(inits)+1||shift.index<=0)
		{
			stop("ERROR:shift.index out of range")
		}
		if(shift.index==1){
			inits<-c(0,inits)
		} else {
			if(shift.index==length(inits)+1){
				inits<-c(inits,0)
			} else {
				#insert in the middle
				inits<-c(inits[1:(shift.index-1)],0, inits[-(1:(shift.index-1))])
			}
		}
		
		#now check the data integrity
		if(length(inits)!=length(num.std))
		{
			stop("ERROR: number of inits doesn't aggree with num.std, please check")
		}
		if(length(x)!=sum(num.std))
		{
			stop("ERROR: number of standards doesn't equal to num.std specified")
		}
		#now seems good the data input
		inits.expand<-rep(inits, num.std);
		xinput<-x+inits.expand
		
		xmid<-pars[1]
		scal<-pars[2]
		g<-1
		if(missing(a))
		{
			a<-min(y)
		}
		if(missing(d))
		{
			d<-max(y)
		}
		
		#pred<-a+(d-a)/(1+exp((xmid-xinput)/scal))^g
		pred<-f5pl(c(a,d,xmid, scal,g),xinput);
		#m<-weight.matrix(pred,model.weight,order)
		(y-pred)  #/m #sqrt(abs(pred))#^2#*(pred)
	}

	##we are looking for a solution to overcome the issue where we run out of limit
## when calculate exp((xmid-x)/scal)
## the solution is like this,
## take a exp(log((d-a)/(1+exp((xmid-x)/scal))^g))
##   rearrange this  exp(log(d-a)-log((1+exp(xmid-x)/scal)^g))
##					exp(log(d-a)-g*log(1+exp(xmid-x)/scal))
##		there are one case we need to consider, when exp is too big or exp is too small (both relative to 1).
##		where exp() is too big or (xmid-x)>709, we approximate to get log(1+exp(xmid-x)/scal)~=(xmid-x)/scal 
##		the whole thing becomes exp(log(d-a)-g*(xmid-x)/scal)
	#the five parameter logistic function
	#take in the parameter 5 parameters and get the function values
	#NOTE: assuming x has been logged 
	#'@title The five-parameter logistic function
	#'@description read in the paramters and the independent variable value(s), and then return the 
	#'		5pl function value(s). For the 4pl model, set g to be 1.
	#'@details The function has the following form \cr
	#'     f(x)=a+(d-a)/((1+exp((xmid-x)/scal))^g)
	#'@param pars numeric the parameters of the 5pl (or 4pl). It has the following content: [a, d, xmid, scal, g].
	#'@param x numeric the log-transformed x value(s).
	#'@return the 5pl function value(s). 
	#'@export
	f5pl<-function(pars,x)
	{
		if(length(pars)!=5){
			stop("pars are not correctly specified!")
		}
		a<-pars[1]
		d<-pars[2]
		xmid<-pars[3]
		scal<-pars[4]
		g<-pars[5]
		
		y<-x; #make the holder first
		y[]<- -1
		#split the array into to parts
		regular<-which((xmid-x)/scal<=709)
		if(length(regular)==0){
			y<-a+exp(log(d-a)-g*(xmid-x)/scal)
		} else {
			y[regular]<-a+(d-a)/(1+exp((xmid-x[regular])/scal))^g
			y[-regular]<-a+exp(log(d-a)-g*(xmid-x[-regular])/scal)	
		}
		return(y)
		#pred<-a+(d-a)/(1+exp((xmid-x)/scal))^g
	}
	
	#'@title Get rid of zeros in a numeric vector 
	#'@description Get rid of zeros in a numeric vector before taking 
	#' the logarithm of them. We basically replace the "zeros" with a negligible small 
	#'	 value in order to avoid NaN upon the log-transformation.
	#'@param x numeric values as input 
	#'@param fac numeric a factor of scale in order to get a "small" value to 
    #' replace zeros 
    #'@return a vector of value with zeros replaced.	
	#'@export
	#Note: output 
	avoidZero<-function(x, fac=10)
	{
		x_min<-0;
		if(sum(x==0)>0)
		{
			x_min<-unique(x)[order(unique(x))[2]]/fac
		}
		return(x+x_min)
	}
	
    #'@title The inverse of the 5-parameter logistic function 
	#'@description The inverse function of the 5pl. Set the value of g to be 1, if the 4pl
	#'	is of interest.
	#'@param pars the parameters of the function. It has the following content: [a, d, xmid, scale, g].
	#'@param y the value to be reverse calculated.
	#'@return the value of the reverse function 
	#'@export
	# #'@seealso [f5pl()]
	# Note: output is exponentiated, but in log scale.
	
	inv.f5pl<-function(pars,y)
	{
		if(length(pars)!=5){
			stop("pars are not correctly specified!")
		}
		a<-pars[1]
		d<-pars[2]
		xmid<-pars[3]
		scal<-pars[4]
		g<-pars[5]
		
		x<-y; #make the holder first
		x[]<- -1
		#split the array into to parts
		#regular<-which((xmid-x)/scal<=709)
		#if(length(regular)==0){
		#	y<-a+exp(log(d-a)-g*(xmid-x)/scal)
		#} else {
		#	y[regular]<-a+(d-a)/(1+exp((xmid-x[regular])/scal))^g
		#	y[-regular]<-a+exp(log(d-a)-g*(xmid-x[-regular])/scal)	
		#}
		
		x[which(y<=a)]=0;
		x[which(y>=d)]=Inf;
		
		#((d-a)/(y-a))^(1/g)
		regular.ind<-which(y>a&y<d);
		yy<-y[regular.ind];
		finite.ind<-which(is.finite(((d-a)/(yy-a))^(1/g)))
		xx<-c();
		if(length(finite.ind)==0)
		{
			xx<-xmid-scal*(1/g)*(log(d-a)-log(yy-a))
		} else {
			xx[finite.ind]<-xmid-scal*log(((d-a)/(yy-a))^(1/g)-1)
			xx[-finite.ind]<-xmid-scal*(1/g)*(log(d-a)-log(yy-a))
		}
		
		xx<-exp(xx);
		x[regular.ind]<-xx;
		return(x)
	}