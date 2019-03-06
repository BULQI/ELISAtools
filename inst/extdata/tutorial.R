#R code to run 5-parameter logistic regression on ELISA data

#load the library
library(ELISAtools)

###
#get file folder
dir_file<-system.file("extdata", package="ELISAtools")

#setwd(dir_file)
batches<-loadData(file.path(dir_file,"design.txt"))

#******IMPORTANT***********
#now set the working directory to somewhere you have permission to write
#*************

#now add
reportHtml(batches,file.dir=tempdir());
#make a guess for the parameters, the other two parameters a and d 
#will be estimated based on data.
model<-"5pl"
pars<-c(7.2,0.5, 0.015) #5pl inits
names(pars)<-c("xmid", "scal", "g")

#model<-"4pl"
#pars<-c(7.2,0.9) #4pl inits
#names(pars)<-c("xmid", "scal")

#do fitting. model will be written into data set.
batches<-runFit(pars=pars,  batches=batches, refBatch.ID=1, model=model  )
#batches<-runFit(pars=pars,  batches=batches, refBatch.ID=1, model="4pl"  )
#now call to do predications based on the model.
batches<-predictAll(batches);

#reporting.
reportHtml(batches, file.name="report_ana",file.dir=tempdir())

#now saving the combine data.
saveDB(batches, file.path(tempdir(),"elisa_tool1.rds"));
batches.old<-loadDB(file.path(tempdir(),"elisa_tool1.rds"));

#now suppose want to join/combine the two batches, old and new 
batches.com<-combineData(batches.old, batches);
reportHtml(batches.com, file.name="report_com",file.dir=tempdir());

batches.com<-runFit(pars=pars,  batches=batches.com, refBatch.ID=1 ,model=model )

#now call to do predications based on the model.
batches.com<-predictAll(batches.com);
reportHtml(batches.com,file.name="report_com_ana", file.dir=tempdir());