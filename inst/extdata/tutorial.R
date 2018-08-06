#R code to run 5-parameter logistic regression on ELISA data

#load the library
library(ELISAtools)

###
#get file folder
dir_file<-system.file("extdata", package="ELISAtools")
setwd(dir_file)
batches<-loadData(file.path(dir_file,"design.txt"))

#make a guess for the parameters, the other two parameters a and d 
#will be estimated based on data.
pars<-c(7.2,0.05, 0.015)
names(pars)<-c("xmid", "scal", "g")

#do fitting. model will be written into data set.
batches<-runFit(pars=pars,  batches=batches, refBatch.ID=1  )

#now call to do predications based on the model.
batches<-predictAll(batches);

#reporting.
reportHtml(batches)