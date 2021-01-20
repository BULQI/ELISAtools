###this is module to take care of input/output for the ELISAtools project
#---- by Feng 07/2018
###############
####import the stringi to take care of locale, mainly in mac 
#'@import stringi
#'@import utils

#'@include ELISAplate.R

#S3 method to quickly annotate the plate with sample ids
#return a 96 well plate with ids. user need to input the sample Id
#but std always be std1 std2 std3
#'@title S3 method to annotate ELISA plate
#'@description to write annotations for an ELISA plate as
#'    an input to guide the functions to read OD values
#' @details Based on the input to quickly write the annotations
#'  for ELISA plate. The output is in a 96-well format and 
#'	will be used to giude the reading of OD plates. This way
#'	only a nxm dataframe can be used. To write non-regular
#'	annotation, you have to do it mannually. 
#'
#'	@param sample.id character vector to specify the names/ids of
#'		the samples on the plate. Note, standard/calibration sample
#'		ids/names is fixed to be "s1","s2", etc, which are specified
#'		by the software and users don't need to privide.
#'	@param sample.prefix characters will be added to the beginning of sample names 
#'	@param sample.suffix characters will be added to the end of sample names
#'	@param num.sample numeric number of samples to write
#'	@param num.std numeric number of standards
#'	@param byRow.sample boolean indicate whether to write sample names 
#'		horizontally by row (TURE) or vertically by column (FALSE)
#'	@param byRow.replicates boolean indicate whether to write sample replicates 
#'		horizontally by row (TURE) or vertically by column (FALSE)
#'	@param replicates.sample numeric number of replicates for each sample
#'	@param replicates.std numeric number of replicates for each standards
#'	@param rows numeric vector to specify which rows to be included in the annotation
#'	@param columns numeric vector to specify which columns 
#'		to be included in the annotation
#'	@param std.first boolean to indicate whether to write standards first or
#'		the samples first.
#' @return a dataframe holding the annotations for the plate.
#'
# #'@examples
# #'
# #'sample.id<-c(1:24)
# #'sample.prefix<-"sam"
# #'sample.suffix<-"_d"
# #'num.sample<-length(sample.id)
# # 'replicates.sample<-2
# # 'replicates.std<-2
# #'num.std<-6
# #'byRow.sample=FALSE;
# #'rows<-c(3:8)
# #'cols<-c(3:12)
# #'ann<-annotate.plate(sample.id=sample.id, sample.prefix=sample.prefix, sample.suffix=sample.suffix,
# #'		num.sample=num.sample, byRow.sample=F, byRow.replicates=T, num.std=num.std,
# #'		rows<-rows, cols<-cols,
# #'		replicates.sample=replicates.sample, replicates.std=replicates.std
# #'		)
# #'
# #'write.table(ann, file=file.path(tempdir(),"annote.txt"), sep="\t", row.names=T,
# #'		col.names=TRUE)
# #' @seealso  \code{\link{SensorgramData-class}} \code{\link{plot}} \code{\link{SaveSPRData}}
# # @export
annotate.plate<-function (sample.id, sample.prefix, sample.suffix,
		num.sample,num.std=8,
		byRow.sample=TRUE, 
		byRow.replicates=TRUE,
		replicates.sample=3, replicates.std=3, 
		rows, columns,std.first=TRUE
		)
{	
	#check the data integraty
	#first see number matches?
	if(missing(sample.id))
	{
		stop("***Error: missing sample id!!");
	}
	len.sample<-length(sample.id);
	row.size<-12; #how many wells a row has?not a total number of rows each plate has
	column.size<-8
	rows.all<-c(1:column.size);
	columns.all<-c(1:row.size);
	if(missing(rows))
	{
		rows<-rows.all;
	} 
	
	if(missing(columns))
	{
		columns<-columns.all;
	}
	num.well<-length(rows)*length(columns);
	#cat("num.well:",num.well,"\n")
	#cat("num fo sample:",num.sample*replicates.sample+num.std*replicates.std,"\n");
	#cat("num std:", num.std,"\n");
	if(num.well!=(num.sample*replicates.sample+num.std*replicates.std))
	{
		stop("Error:sample number and standard number do not equal to \n\tnumber of wells on the plate, please check!!")
	}
	#now populate the plate data frame
	#row.index<-rows
	row.name<-c("A","B","C", "D", "E","F", "G","H")
	row.name<-row.name[rows]
	#column.index<-c(1:row.size)
	#column.index<-column.index[-exclude.column]
	col.name<-c("1","2","3", "4", "5","6", "7","8", "9","10","11","12")
	col.name<-col.name[columns]
	
	#prepare the sample annotations
	
	if(!missing(sample.prefix))
	{
		sample.id<-paste0(sample.prefix, sample.id)
	}
	if(!missing(sample.suffix))
	{
		sample.id<-paste0( sample.id, sample.suffix)
	}
	#add in std 
	std.id<-paste0("s",c(1:num.std))
	#byRow.replicates<-match.arg(byRow.replicates)
	#if(byRow.replicates){
	#	#horizontally
	#	std.id<-rep(std.id,rep(replicates.std,num.std))
	#	sample.id<-rep(sample.id,rep(replicates.sample, length(sample.id)))
	#}
	#sample.ids<-c(std.id,sample.id)
	#std.first<-match.arg(std.first)
	#if(!std.first){
	#	sample.ids<-c(sample.id, std.id)
	#}
	plate<-NULL
	#byRow.sample<-match.arg(byRow.sample)
	if(byRow.sample){
		if(byRow.replicates){
			#make the sample id & std id ready by repeating and concatenating them
			std.id<-rep(std.id,rep(replicates.std,num.std))
			sample.id<-rep(sample.id,rep(replicates.sample, length(sample.id)))
			sample.ids<-c(std.id,sample.id)
			#std.first<-match.arg(std.first)
			if(!std.first){
				sample.ids<-c(sample.id, std.id)
			}
			#building the data.frame
			for(i in 1:length(rows))
			{
				plate<-rbind(plate,sample.ids[c(1:length(columns))+(i-1)*length(columns)])
			}
		} else { #in this case, we also need to replicate elements first, but vertically. 
			##luckily in this case, the repeats and well number are identical.
			
			#this case is bit difficult, we need to check more
			if((replicates.std!=replicates.sample)&&(num.std%%length(columns)!=0)&&(num.sample%%length(columns)!=0))
			{
				stop("ERROR: incompatible format: 1)non-equal replicate number between standard and sample;2)std not cover the whole row.");
			}
			#index<-1;
			rep.num<-0;
			sample.ids<-c(std.id,sample.id)
			if(!std.first){
				sample.ids<-c(sample.id, std.id)
			}
			if(replicates.std==replicates.sample)  ##replicates numbers are identical, we just go ahead to replicate rows
			{
				rep.num<-rep(replicates.sample, length(sample.ids)/length(columns))
			} else { #in this case, replicates between std and sample are not equal. tricky, we need to do more to figure out	
					#good things though, num.std must be a multiple of the row numbers
				#figure out rep number for each row
				rep.num<-c(rep(replicates.std,num.std/length(columns)),
						rep(replicates.sample, num.sample/length(columns)));
				if(!std.first)
				{
					rep.num<-c(rep(replicates.sample, num.sample/length(columns)),
						rep(replicates.std,num.std/length(columns)));
				}
			} #done.
						
			for(i in 1:length(rep.num))
			{
				for(j in c(1:rep.num[i]))
				{
					plate<-rbind(plate,sample.ids[c(1:length(columns))+(i-1)*length(columns)])
				}
			}
		}
	} else {  ###vertically fill the sample
		if(byRow.replicates){ ##this is a harder one, we need first to check for consistency, this is redundant 
			#this case is bit difficult, we need to check more
			if((replicates.std!=replicates.sample)&&(num.std%%length(rows)!=0)&&(num.sample%%length(rows)!=0))
			{
				stop("ERROR: incompatible format: 1)non-equal replicate number between standard and sample;2)std not cover the whole column.");
			}
			rep.num<-0;
			sample.ids<-c(std.id,sample.id)
			if(!std.first){
				sample.ids<-c(sample.id, std.id)
			}
			if(replicates.std==replicates.sample)
			{
				rep.num<-rep(replicates.sample, length(sample.ids)/length(rows));
			} else { #again this is hard, but the good thing is that we have the number of std is a multiply of column length
			#figure out rep number for each row
				rep.num<-c(rep(replicates.std,num.std/length(rows)),
						rep(replicates.sample, num.sample/length(rows)));
				if(!std.first)
				{
					rep.num<-c(rep(replicates.sample, num.sample/length(rows)),
						rep(replicates.std,num.std/length(rows)));
				}
			}
			for(i in 1:length(rep.num))
			{
				for(j in c(1:rep.num[i]))
				{
					plate<-cbind(plate,sample.ids[c(1:length(rows))+(i-1)*length(rows)]);
				}
			}
		} else { #this should be easy???
			#make the sample id & std id ready by repeating and concatenating them
			std.id<-rep(std.id,rep(replicates.std,num.std))
			sample.id<-rep(sample.id,rep(replicates.sample, length(sample.id)))
			sample.ids<-c(std.id,sample.id)
			#std.first<-match.arg(std.first)
			if(!std.first){
				sample.ids<-c(sample.id, std.id)
			}
			#building the data.frame
			for(i in 1:length(columns))
			{
				plate<-cbind(plate,sample.ids[c(1:length(rows))+(i-1)*length(rows)])
			}
		}
		
		
	} #end of filling sample byRow or not
	
	rownames(plate)<-row.name;
	colnames(plate)<-col.name;
	
	return(plate);
}

#reading the annotation for the plate
#this will give us the information about sample identities
#'@title Read the annotation of single ELISA plate
#'@description Parse the annotations for one single ELISA plate from
#'    a section of a file and output the annotations for standard and unknown
#'		separately.
#' @details The annotation file may contain annotations for more than
#'		one plate. Each plate is marked by "Plate: plate 1..." and 
#'		"~End". This function is fed in with the content for each section
#'		and we do actually parsing in here. Store the annotations into 
#'		data frame. It also parse the standard concentration and 
#'		include this information in the data frame.
#'		For each section,we expect
#'		the following format\cr
#'	\tabular{llllllll}{
#'    \tab	1\tab	2\tab		3\tab		4\tab 			...\cr	
#'   C\tab	s1\tab	s1\tab		sample1\tab		sample1\tab	...\cr
#'   D\tab	s2\tab	s2\tab		sample2\tab		sample\tab  ...\cr	
#'	...\tab	...\tab	...\tab		...\tab		...\tab		...\cr
#'	}
#'		In addition, the row name and column names indicate the
#'		the plate row and column indices.
#'		As input, the stardard and unknown are returned separately in two 
#'		tables.
#'	@param annotation characters to specify the path and name of the annotation file
#'	@param std.conc data.frame containing standard concentration data. 
#'		Only first two columns are used with first one to be the standard IDs
#'		and second the concentrations.
# #'	@param file.dir file path to the annotatoin file.
#' @return a list of data.frames holding the annotations for the plate.
#'
#'@examples
#'#get example annotation file path from the system folder
#' fileName<-system.file("extdata", "annote_single.txt", package="ELISAtools")
#'#prepare the standard concentration file.
#' std.conc<-data.frame(id=c("s1","s2","s3","s4","s5","s6"), conc=c(1:6))

#'#read the data as a data frame.
#' ann<-read.table(fileName, header=TRUE,  sep="\t", stringsAsFactors=FALSE)
#'
#'#call to do the reading.
# #' annotation<-read.annotation(ann,  std.conc)
# #' @seealso  \code{\link{elisa_batch-class}} \code{\link{loadData}} 
# #' @export
read.annotation<-function(annotation,  std.conc)
{
	#check the data integraty 
	if(missing(annotation))
	{
		stop("ERROR: no annotation input specified");
	}
	if(class(annotation)!="data.frame")
	{
		stop("ERROR:  annotation input is not in a correct format");
	}
#	if(!file.exists(annotation))
#	{
#		stop("ERROR: annotation does not exist, please check.");
#	}
#	if(missing(file.dir))
#	{
#		file.dir<-file.path("./")
#	}
#	ann<-read.table(file.path(file.dir, annotation), header=T, sep="\t", stringsAsFactors=F);
	ann<-annotation
	rows<-rownames(ann);
	#rows<-as.numeric(rows)
	rows.ind<-c(1:8)
	names(rows.ind)<-c("A","B","C","D","E","F","G","H")
	if(!all(is.element(rows, names(rows.ind))))
	{
		stop("ERROR:the row index is not correct set, please check!")
	}
	rows<-rows.ind[rows]
	cols<-colnames(ann)
	cols<-sub("^[A-Z]+", "", cols,ignore.case=T)
	cols<-as.numeric(cols)
	if(!all(cols<=12 & cols>=1))
	{
		stop("ERROR:the col index is not correct set, please check!")
	}

	#if(missing(rows))
	#{
	#	rows<-c(1:8)
	#}
	#if(missing(cols))
	#{
	#	cols<-c(1:12)
	#}
	
	if(dim(ann)[1]>8||dim(ann)[2]>12)
	{
		stop("ERROR: annotation data frame is not in a correct format.");
	}
	
	if(missing(std.conc))
	{
		stop("ERROR:missing standard concentration data, please specify")
	}
	if(class(std.conc)!="data.frame")
	{
		stop("ERROR:std.conc should be a data frame");
	}
	if(dim(std.conc)[2]!=2)
	{
		warning("There are more columns in standard concentration data table, only first two are used",
			immediate.=TRUE);
	}
		
	#start parsing
	#sample.ids<-(as.vector(ann))
	#std.idx<-grep("^s[1-9]+$",sample.ids)
	#std.ids<-(sample.ids[std.idx])
	#sample.ids<-(sample.ids[-std.idx])
	
	annotations.std<-NULL #data.frame(ID=std.ids,row=0,col=0,conc=0)
	annotations.unknown<-NULL #data.frame(ID=std.ids,row=0,col=0)
	
	#parse and read row&col
	#std first
	count.std<-1
	count.unknown<-1
	for(i in 1:dim(ann)[1]) #row
	{
		row<-rows[i];
		#cat("i:",i,";\n")
		for(j in 1:dim(ann)[2]) #col
		{
			col<-cols[j];
			#cat("\tj:",j,";")
			#"NA" 
			if(is.na(ann[i,j]))
			{
				next;
			}
			#white space
			if(length(grep("^\\s*$",ann[i,j]))==1)
			{
				next;
			}
			if(length(grep("^masked|mask$",ann[i,j],ignore.case=T))==1)
			{
				next;
			}
			if(length(grep("^s[1-9]+$",ann[i, j]))==1)
			{
				#cat("\t in std\n");
				if(count.std==1)
				{
					annotations.std<-data.frame(ID=ann[i,j],row=row,col=col,conc=0, stringsAsFactors=F,row.names=NULL)
				} else {
					annotations.std[count.std,1]<-ann[i,j];#<-rbind(annotations.std,data.frame(ID=ann[i,j], row;
					annotations.std[count.std,2]<-row;
					annotations.std[count.std,3]<-col;
				}
				#look up the concentration
				ind<-which(std.conc[,1]==ann[i,j])
				if(length(ind)==0)
				{
					stop("ERROR:can not find the concentration for the std, please check!!")
				}
				if(length(ind)>1)
				{
					warning("more than one concentration specified for std sample, only the first one is used!!") 
				}
				
				annotations.std[count.std,4]<-std.conc[ind[1],2]
				count.std<-count.std+1;
			} else { #sample
				#cat("\tNOT\n")
				if(count.unknown ==1)
				{
					annotations.unknown<-data.frame(ID=ann[i,j],row=row, col=col, stringsAsFactors=F,row.names=NULL)
				} else {
					annotations.unknown[count.unknown,1]<-ann[i,j];
					annotations.unknown[count.unknown,2]<-row;
					annotations.unknown[count.unknown,3]<-col;
				}
					
				count.unknown<-count.unknown+1;
			}
		} 
	}
	if(is.null(annotations.unknown)){
		annotations.unknown<-data.frame();
	}else{
		
		class(annotations.unknown$col)<-"integer"; 
	}
	class(annotations.std$col)<-"integer";
	list(standards=annotations.std, unknowns=annotations.unknown)
}
#
#reading the annotations for plates
#this will give us the information about sample identities
#'@title Read the annotations of plates 
#'@description Parse annotations for multiple ELISA plates from
#'   	files, one annotation file and one standard concentration file,
#'		and output the annotations for standard and unknown
#'		separately.
#' @details The annotation file may contain annotations for more than
#'		one plate. Each plate is marked by "Plate: plate 1..." and 
#'		"~End". This function parses each section in both annotation
#'		file and standard concentration file. Then passes the section
#'		on to do the parsing.
#'		For each section,we expect
#'		the following format\cr
#'	\tabular{llllllll}{
#'    \tab	1\tab	2\tab		3\tab		4\tab 			...\cr	
#'   C\tab	s1\tab	s1\tab		sample1\tab		sample1\tab	...\cr
#'   D\tab	s2\tab	s2\tab		sample2\tab		sample\tab  ...\cr	
#'	...\tab	...\tab	...\tab		...\tab		...\tab		...\cr
#'	}
#'	@param annotation characters to specify the path and name of the annotation file
#'	@param std.conc characters to specify the standard concentration file. 
#'		
#'	@param dir.annotation characters specifying the file to the annotatoin file.
#'	@param dir.stdConc characters specifying the path to the annotatoin file.
#'	@param num.plate numeric indicating the number of plates in the annotation
#'			files.
#' @return a list of annotations for elisa plates.
#'
#'@examples
#'#get example annotation file path from the system folder
#' ann<-system.file("extdata", "annote.txt", package="ELISAtools")
#' std.conc<-system.file("extdata", "stdConc.txt", package="ELISAtools")
#'
#'#read them in and there are 2 plates.
#'	read.annotations(annotation=ann,  std.conc=std.conc, num.plate=2)
# #' @seealso  \code{\link{elisa_run-class}} \code{\link{loadData}} 
#' @export
read.annotations<-function(annotation,  std.conc, dir.annotation, dir.stdConc,num.plate=1)
{
	#check the data integraty 
	if(missing(annotation))
	{
		stop("ERROR: no annotation input specified");
	}
	if(!missing(dir.annotation))
	{
		annotation<-file.path(dir.annotation, annotation)
	}
	
	#annotation<-file.path(dir.annotation, annotation)
	if(!file.exists(annotation))
	{
		stop("ERROR:  annotation input does not exists");
	}
	if(missing(std.conc))
	{
		stop("ERROR: no std conc input specified");
	}
	
	if(!missing(dir.stdConc))
	{
		std.conc<-file.path(dir.stdConc, std.conc);
	}
	if(!file.exists(std.conc))
	{
		cat("std.conc:",std.conc,"\n");
		stop("ERROR: annotation does not exist, please check.");
	}
	
	#read the file first
	con.ann<-file(annotation,"r");
	ann.raw<-readLines(con.ann, skipNul=T); #read all
	close(con.ann);
	con.stdConc<-file(std.conc, "r")
	stdConc.raw<-readLines(con.stdConc, skipNul=T); #read all
	close(con.stdConc);

	#first check to see whether there are enough elements
	anns<-list();
	for( i in 1:num.plate)
	{
		#find the anchors or markers for each section/plate
		ind.start<-grep("^Plate:",ann.raw)
		if(length(ind.start)<1)
		{
			warning("not enough plates found in annotation file, please check");
			cat("not enough plates found in annotation file, please check\n");
			
			break;
			
		}
		ind.start<-ind.start[1]+1;
		ind.end<-grep("^~End",ann.raw);
		if(length(ind.end)<1)
		{
			warning("not enough plates found in annotation file, please check");
			ind.end<-length(ann.raw);
		}
		ind.end<-ind.end[1]-1
		ann<-ann.raw[c(ind.start:ind.end)];
		ann.raw<-ann.raw[c((ind.end+1+1):length(ann.raw))];
		
		#doing standard
		ind.start<-grep("^Plate:", stdConc.raw);
		if(length(ind.start)<1)
		{
			warning("not enough plates found in standard conc file, please check");
			break;
		}
		ind.start<-ind.start[1]+1;
		ind.end<-grep("^~End",stdConc.raw);
		if(length(ind.end)<1)
		{
			warning("not enough plates found in standard conc file, please check");
			ind.end<-length(stdConc.raw);
		}
		ind.end<-ind.end[1]-1
		sConc<-stdConc.raw[c(ind.start:ind.end)];
		stdConc.raw<-stdConc.raw[c((ind.end+1+1):length(stdConc.raw))];
		#now we need to read the second into data frame 
		ann.df<-read.table(text=ann, sep="\t", stringsAsFactors=F, header=T, row.names=1);
		sConc.df<-read.table(text=sConc, sep="\t", stringsAsFactors=F, header=T);
		
		#feed to do parsing
		anns[[i]]<-read.annotation(annotation=ann.df, std.conc=sConc.df);
	}
	return(anns)
}

#read each individual file of ODs 
#reading the annotation for the plate
#this will give us the information about sample identities
#'@title Read the single ELISA OD plate
#'@description Read the individual ELISA plate to parse the ODs.
#'     
#'@details The input is a text file imported from the sdf file.  
#'	We only read the first section with both the OD and blank file.
#'	The OD data are read in according to the annotation file.
#'	
#'	@param ODs characters containing data of ODs for one plate
#'	@param annotation list of data containing annotations of the plate
#'	@param batchID characters specifying the batchID read from the design file
#'	@param expID characters specifying the expID or plateID read from the design file 
#' @return an object of elisa_plate holding data and annotations for 
#'		a single plate.
#'
# #'@examples
# #' #setwd( "E:\\feng\\LAB\\hg\\ELISA\\ELISAtools\\dev")
# #' std.conc<-data.frame(id=c("s1","s2","s3","s4","s5","s6"), conc=c(1:6))
# #' annotation<-"annote.txt"
##'
# #'	read.annotation(annotation,  std.conc)
# #' @seealso  \code{\link{elisa_batch-class}} \code{\link{loadData}} \code{\link{read.plates}}
# #' @export
read.plate<-function(ODs, annotation, batchID, expID)
{
	if(missing(ODs))
	{
		stop("Error: please specify the files of OD data")
	}
	if(class(ODs)!="character")
	{
		stop("ERROR: the input OD is not in correct format!!")
	}
	if(missing(annotation))
	{
		stop("ERROR: please specify the annoation data for the plate");
	}
	if(class(annotation)!="list")
	{
		stop("ERROR: the input annotation is not in correct format")
	}
	if(is.null(annotation$standards)||class(annotation$standards)!="data.frame")
	{
		stop("ERROR: the input annotation is not in correct format")
	}
	if(!is.null(annotation$unknowns)&&class(annotation$unknowns)!="data.frame")
	{#could be null, meaning no unknowns
		stop("ERROR: the input annotation is not in correct format")
	}
	
	###now read in the plate 
	#con<-file(fileName,"r");
	#OD.raw<-readLines(con, skipNul=T,n=100); #read in 100 lines, should be way more than enough
	#close(con);
	OD.desc<-NULL
	OD.header<-NULL
	OD.plate<-NULL
	OD.blank<-NULL
	OD.readOnce<-FALSE
	for(i in 1:length(ODs))
	{
		#cat("i:",i,"\n")
		tempStr<-ODs[i];#trimws(ODs[i], which="both");
		if(length(grep("^[ \t\r\n]*~End",tempStr))>0)
		{#we are done.
			#cat("\tbreak");
			break;
		}
		if(length(grep("^Plate", tempStr))>0)
		{
			OD.desc<-tempStr;
			next;
		}
		if(length(grep("Temperature", tempStr))>0)
		{
			OD.header<-tempStr;
			next;
		}
		#if(length(grep("^Plate", tempStr))>0)
		#{
		#	OD.desc<-tempStr;
		#	next;
		#}
		if(nchar(trimws(tempStr,which="both"))>0)
		{
			if(!OD.readOnce)
			{
				OD.plate<-tempStr;
				OD.readOnce<-TRUE;
			} else {
				OD.blank<-tempStr;
			}
			next;
		}
	}
	#now we have everything read in.
	#start parsing.
	#leave the description not parsed
	#header
	OD.header<-strsplit(OD.header, "\t")[[1]]
	OD.header<-OD.header[-c(1:2)];
	
	OD.plate<-strsplit(OD.plate, "\t")[[1]]
	temperature<-as.numeric(OD.plate[2])
	OD.plate<-as.numeric(OD.plate[-c(1:2)]);
	if(is.null(OD.blank))
	{
		OD.blank<-OD.plate;
		OD.blank[]<-0;
	}	else {
		OD.blank<-strsplit(OD.blank, "\t")[[1]]
		#temperature<-as.numeric(OD.blank[2])
		OD.blank<-as.numeric(OD.blank[-c(1:2)]);
	}
	
	#now we are done parsing, put into data frame and get ready to do output
	eplate<-elisa_plate();
	#cat("OD.header length:",length(OD.header), "\n");
	eplate<-load.ODs(eplate,plate.header= OD.header, plate.data=OD.plate, 
		plate.blank=OD.blank, annotation=annotation);
	#cat("\tdone.......\n");
	eplate@batchID<-batchID;
	eplate@expID<-expID;
	eplate@desc<-OD.desc;
	eplate@range.ODs<-c(min(eplate@data.std$OD),max(eplate@data.std$OD));
	return(eplate)
}

#'@title Read the ELISA OD files
#'@description Read the ELISA OD file to parse the ODs.
#'     
#'@details The input is a text file imported from the sdf file.  
#'	The file may contain multiple plates of OD. We will parse 
#'	each file section and then read them according to the 
#'	annotation to load the data. We assume for each file the data
#'	are for the same batch and experiment. If otherwise, please 
#'	split the file into different ones.
#'	
#'	@param fileName characters containing file name of OD data
#'	@param annotations list of data containing annotations of the plates
#'	@param batchID characters specify the batchID read from the design file
#'	@param expID characters specify the expID or plateID read from the design file
#'	@param num.plate numeric number of OD plates in the OD file. 
#'	@param date characters the date running the ELISA exps.
#' @return an object of elisa_run holding data and annotations for 
#'		one or multiple plates.
#'
#'@examples
#'#get example annotation file path from the system folder
#' ann<-system.file("extdata", "annote.txt", package="ELISAtools")
#' std.conc<-system.file("extdata", "stdConc.txt", package="ELISAtools")
#'
#'#read them in and there are 2 plates.
#'	annotations<-read.annotations(annotation=ann,  std.conc=std.conc, num.plate=2)
#'
#'#now start reading the OD plate file
#' fileName <-system.file("extdata", "Assay_3_and_4.txt", package="ELISAtools")
#'	plates<-read.plates(fileName, annotations=annotations, num.plate=2, batchID="b1", expID="e1")
# #' @seealso  \code{\link{elisa_batch-class}} \code{\link{loadData}} \code{\link{read.plates}}
#'
#'@export
read.plates<-function(fileName, annotations, num.plate=1, batchID, expID ,date=NA_character_)
{
	if(missing(fileName))
	{
		stop("Error: please specify the files of OD data")
	}
	if(!file.exists(fileName))
	{
		stop("ERROR: can not find the specified file!!")
	}
	if(missing(annotations))
	{
		stop("ERROR: please specify the annoation data for the plate");
	}
	if(class(annotations)!="list")
	{
		stop("ERROR: the input annotation is not in correct format")
	}
	if(missing(batchID))
	{
		stop("ERROR: please specify the batch ID");
	}
		if(missing(expID))
	{
		stop("ERROR: please specify the exp ID");
	}
	#if(is.null(annotation$standards)||class(annotation$standards)!="data.frame")
	#{
	#	stop("ERROR: the input annotation is not in correct format")
	#}
	#if(is.null(annotation$unknowns)||class(annotation$unknowns)!="data.frame")
	#{
	#	stop("ERROR: the input annotation is not in correct format")
	#}
	
	###now read in the plate 
	con<-file(fileName,"r");
	OD.raw<-readLines(con, skipNul=T); #read in 100 lines, should be way more than enough
	close(con);
	OD.raw<-stri_conv(OD.raw, to="UTF-8");
	if(length(annotations)<num.plate)
	{
		stop("ERROR: the num of plate specified is less than the elements of annotatoins, please check");
	}

	eplates<-elisa_run();
	eplates@num.plates<-num.plate
	eplates@batchID<-batchID
	eplates@desc<-expID;
	eplates@date<-date;
	#OD.desc<-NULL
	#OD.header<-NULL
	#OD.plate<-NULL
	#OD.blank<-NULL
	#OD.readOnce<-FALSE
	range.ODs.min<- 100;
	range.ODs.max<- 0
	for(i in 1:num.plate)
	{
		#cat("&&&&&&&&&plate:", i, "\n")
		#find the anchors or markers for each section/plate
		ind.start<-grep("^Plate:",OD.raw)
		if(length(ind.start)<1)
		{
			warning("not enough plates found in OD data file (no starting point), please check");
			#stop("not enough plates found in OD data file (no starting point), please check");
			
			break;	
		}
		ind.start<-ind.start[1];
		ind.end<-grep("^[ \t]*~End[ \t\r\n]*$",OD.raw, ignore.case=T);
		if(length(ind.end)<1)
		{
			warning("not enough plates found in OD data file (no ending point), please check");
			#stop("not enough plates found in OD data file (no ending point), please check");
			ind.end<-length(OD.raw);
		}
		ind.end<-ind.end[ind.end>ind.start][1]
		OD<-OD.raw[c(ind.start:ind.end)];
		OD.raw<-OD.raw[c((ind.end+1):length(OD.raw))];
		
		#now feed to read plate
		eplates@plates[[i]]<-read.plate(OD, annotations[[i]], batchID,expID);
		if(eplates@plates[[i]]@range.ODs[1]<range.ODs.min)
		{
			range.ODs.min<-eplates@plates[[i]]@range.ODs[1];
		}
		if(eplates@plates[[i]]@range.ODs[2]>range.ODs.max)
		{
			range.ODs.max<-eplates@plates[[i]]@range.ODs[2]
		}
		
	}
	
	#now we are done parsing, put into data frame and get ready to do output
	#eplate<-elisa_plate();
	#eplate<-load.ODs(eplate,plate.header= OD.header, plate.data=OD.plate, plate.blank=OD.blank, annotation=ann);
	#eplate@batchID<-batchID;
	#eplate@expID<-expID;
	#eplate@desc<-OD.desc;
	eplates@range.ODs<-c(range.ODs.min, range.ODs.max)
	return(eplates)
}


#####read in 
##for annotation, we assume there is header indicating the row index,and rownames for row index correctly!!.
#'@title Read data according to the design file
#'@description Read the design file and then load the 
#'	data according to the information in the design file.
#'     
#'@details The design file contains all the information necessary
#'	to read data. It has the following format\cr
#'	\tabular{lllllllll}{
#'	ExpID\tab	FileName\tab	Batch\tab	Num_Plate\tab	Date\tab	AnnotationFile\tab	Std_Conc\tab	Dir_Annotation\tab	Dir_StdConc\cr
#'	Exp1\tab	file1.txt\tab	Batch1\tab	1\tab	9/18/2009\tab	annote.txt\tab	stdConc.txt\tab \tab \cr		
#'	Exp2\tab	file2.txt\tab	Batch2\tab	2\tab	9/18/2009\tab	annote.txt\tab	stdConc.txt\tab \tab \cr		
#'	}
#'	The return data is a list of batches (\code{\link{elisa_batch-class}}), 
#'	which are made of one or many elisa runs(\code{\link{elisa_run-class}}) 
#'	. The
#'		run could contain one or many elisa plates (\code{\link{elisa_plate-class}})
#'		with data or annotation of each plate.\cr
#'  \tabular{lllllll}{
#'	list\tab	|\tab	 \tab	 \tab	 \tab	 \tab	 \cr
#'	 \tab	|\tab	--batch1\tab	|\tab	\tab	\tab	\cr		
#'	 \tab	|\tab	 \tab	|\tab	--run1\tab	|\tab	\cr		
#'	 \tab	|\tab	 \tab	|\tab	 \tab	|\tab	--plate1\cr
#'	 \tab	|\tab	 \tab	|\tab	 \tab	|\tab	--plate2\cr
#'	 \tab	|\tab	--batch2\tab	|\tab	\tab	\tab	\cr
#'	 \tab	|\tab	 \tab	|\tab	\tab	\tab	\cr				
#'}
#'			
#'	@param  design.file characters to specify the path and the file name of the design file.
#' @return a list of batches holding different runs of elisa, which could contain
#'		one or many elisa_plates with data and annotations for 
#'		each plate.
#'
#'@examples
#' file.dir<-system.file("extdata", package="ELISAtools")
#' loadData(file.path(file.dir,"design.txt"));
#' @seealso  \code{\link{elisa_batch-class}} \code{\link{elisa_plate-class}} \code{\link{elisa_run-class}}
#'
#'@export
loadData<-function(design.file)
{
	#first check the data integrity
	if(missing(design.file))
	{
		stop("ERROR: design file is missing, please specify");
	}
#	if(missing(dir.annotation))
#	{
#		dir.annotation<-file.path(".")
#	}
#	if(missing(dir.ODs))
#	{
#		dir.ODs<-file.path(".")
#	}
#	if(missing(dir.std))
#	{
#		dir.std<-file.path(".")
#	}
	
	dfile<-read.table(design.file, header=T, sep="\t", stringsAsFactors=F);
	
	#start reading each individual plate
	#into different batch
	batch.IDs<-unique(dfile$Batch)
	batches<-list();#vector(mode="list",length=length(batch.IDs))
	dir_design<-normalizePath(dirname(design.file));
	for(i in 1:length(batch.IDs))
	{
		cat("Reading Data for Batch:",i,"--",batch.IDs[i],"\n")
		ind<-which(dfile$Batch==batch.IDs[i])
		ebatch<-elisa_batch();
		ebatch@batchID<-batch.IDs[i];
		ebatch@num.runs<-length(ind);
		range.ODs.min<-1000;
		range.ODs.max<-0;
		for(j in 1:length(ind))
		{
			cat("\tExperiment:",j,"--", dfile[ind[j],]$ExpID,"\n")
			flush.console();
			
			dir_ann<-dir_design;
			
			if(!is.na(dfile[ind[j],]$Dir_Annotation))
			{
				dir_ann<-normalizePath(dfile[ind[j],]$Dir_Annotation);
			}
			if(dirname(dfile[ind[j],]$AnnotationFile)!=".")
			{
				dir_ann<-normalizePath(dirname(dfile[ind[j],]$AnnotationFile));
			}			
			dir_sc<-dir_design
			
			if(!is.na(dfile[ind[j],]$Dir_StdConc))
			{
				dir_sc<-normalizePath(dfile[ind[j],]$Dir_StdConc);
			}
			if(dirname(dfile[ind[j],]$Std_Conc)!=".")
			{
				dir_sc<-normalizePath(dirname(dfile[ind[j],]$Std_Conc));
			}
			#need to read annotation
			annotations<-read.annotations(annotation=dfile[ind[j],]$AnnotationFile, 
					std.conc=dfile[ind[j],]$Std_Conc, num.plate=dfile[ind[j],]$Num_Plate,
					dir.annotation=dir_ann, 
					dir.stdConc=dir_sc);
			    
			ebatch@runs[[j]]<-read.plates(fileName=file.path(dir_design,dfile[ind[j],]$FileName), 
					annotations=annotations, 
					batchID=dfile[ind[j], ]$Batch, expID=dfile[ind[j],]$ExpID,
					num.plate=dfile[ind[j],]$Num_Plate,date=dfile[ind[j],]$Date
					)
			#cat("in here")
			#ebatch@runs[[j]]@date<-
			if(ebatch@runs[[j]]@range.ODs[1]<range.ODs.min)
			{
				range.ODs.min<-ebatch@runs[[j]]@range.ODs[1];
			}
			if(ebatch@runs[[j]]@range.ODs[2]>range.ODs.max)
			{
				range.ODs.max<-ebatch@runs[[j]]@range.ODs[2]
			}
		}
		ebatch@range.ODs<-c(range.ODs.min, range.ODs.max);
		batches[[ batch.IDs[i] ]]<-ebatch;
	}
	cat("Done!!!\n")
	return(batches)
}

#functions to 
#'@title Read the saved elisa_batch data
#'@description Load the serialized elisa_batch data from disk.
#'     
#'@details Here we deserialize elisa_batch data by wrapping the readRds()
#'		function call.
#'		The serialized elisa_batch data are assumed to have been correctly 
#'		analyzed. We will print a summary for what has been read. 
#'@param  db characters to specify the path and file name the elisa data file.
#' @return a list of batches holding different runs of elisa, which could contain
#'		one or many elisa_plates with data and annotations for 
#'		each plate.
#'
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
#'#now call to do predications based on the model.
#'batches<-predictAll(batches);
#'
#'#now saving the data.
#'saveDB(batches, file.path(tempdir(),"elisa_tool1.rds"));
#'
#' loadDB(file.path(tempdir(),"elisa_tool1.rds"));
#' @seealso  \code{\link{elisa_batch-class}} \code{\link{loadData}} \code{\link{saveDB}}
#'
#'@export
loadDB<-function(db)
{
	if(missing(db))
	{
		stop("ERROR:please specify the input");
	}
	if(!file.exists(db))
	{
		if(substr(db,nchar(db)-3, nchar(db))!=".rds")
		{
			db<-paste0(db,".rds")
		}
		if(!file.exists(db))
		{
			stop("ERROR:Can not find the specified ELISA database, Please specify!!")
		}
	}
	cat("  ***loading ELISA data set: ",db,"\n")
	ret<-readRDS(db);
	cat("  ***Success!!\n");
	return(ret)
}

#function to serialize the elisa_batch data 

#functions to 
#'@title Save the elisa_batch data
#'@description Serialize elisa_batch data to disk.
#'     
#'@details We serialize elisa_batch data by wrapping the saveRds()
#'		function call.
#'		The serialized elisa_batch data are assumed to have been correctly 
#'		analyzed. We will print a summary for what has been saved. 
#'@param  db character the file name specifying name of the db.
#'@param  batches list of elisa batch data to be serialized.
#'
#'
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
#'#now call to do predications based on the model.
#'batches<-predictAll(batches);
#'
#'#now saving the data.
#'saveDB(batches, file.path(tempdir(),"elisa_tool1.rds"));
#'
#' @seealso  \code{\link{elisa_batch-class}} \code{\link{loadData}} \code{\link{saveDB}}
#'
#'@export
saveDB<-function(batches, db)
{
	if(missing(batches)||missing(db))
	{
		stop("one or both of the input missing, please check!!") 
	}
	if(substr(db,nchar(db)-3, nchar(db))!=".rds")
	{
		db<-paste0(db,".rds")
	}
	if(file.exists(db))
	{
		cat("the specified database file exists, and will be overwritten");
	}
	cat("  ***saving ELISA data set: ",db,"\n")
	ret<-saveRDS(batches, db);
	cat("  ***success!!\n");
	return(ret);
}	

#functions to 
#'@title Save elisa_batch analysis results
#'@description Save the data analysis results to disk in text format.
#'     
#'@details The results are written to disk in the text format (tab-delimited) and is 
#'	easy to be used for other analysis. 
#'@param  file.name character specifying name of the output file.
#'@param  batches list of elisa batch data to be serialized.
#'
#'
#'@examples
#' #'#R code to run 5-parameter logistic regression on ELISA data
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
#'#now call to do predications based on the model.
#'batches<-predictAll(batches);
#'
#'#now saving the data in text.
#'saveDataText(batches, file.path(tempdir(),"elisa_data.txt"));
#'
# #' @seealso  \code{\link{elisa_batch-class}} \code{\link{loadData}} \code{\link{saveDB}}
#'
#'@export
saveDataText<-function(batches, file.name)
{
	if(missing(batches))
	{
		stop("please specify the input batch data");
	}
	if(missing(file.name))
	{
		stop("please specify the file name for the analysis results");
	}
	if(file.exists(file.name))
	{
		cat("the specified file for saving analysis results exists. It will be overwritten");
	}
	file.conn<-file(file.name);
	open(file.conn,open="w");
	write(c("ELISA tool data analysis results"),file.conn, append=FALSE);
	write(c(paste0("Date:\t",date(),"\r\n")),file.conn, append=TRUE);
	
	for( i in 1:length(batches))
	{
		write(c("==============="),file.conn, append=TRUE)
		
		batch<-batches[[i]];
		#for each batch write the following.........
		write(c(paste0("batch:\t",batch@batchID,"\tS Factor:\t",batch@normFactor)),file.conn, append=TRUE)
		write(c("==============="),file.conn, append=TRUE)
		#now we first rearrange the data into data.frame and the write it out.
		for(j in 1:batch@num.runs)
		{
			#write(c(paste0("R:\t",batch@batchID,"\tS Factor:\t",batch@normFactor)),file.conn
			for(k in 1:batch@runs[[j]]@num.plates)
			{
				write(paste0("RUN_#",j,"\t",batch@runs[[j]]@date, "\tplate_#",k,"\tS Factor:",batch@runs[[j]]@plates[[k]]@normFactor),file.conn, append=TRUE);
				
				suppressWarnings(write.table(batch@runs[[j]]@plates[[k]]@mdata.std,file=file.conn,append=T,sep="\t", row.names = F))
				#start making the data frame for output data
				unknown<-batch@runs[[j]]@plates[[k]]@data.unknown;
				if(is.null(unknown)||dim(unknown)[1]==0){
					next;
				}
				ids<-unique(unknown$ID);
				#determine the number repeats 
				nRep<-max(aggregate(unknown, FUN=length,by=list(unknown$ID))$ID)
				
				dfm<-data.frame();
				if(!is.null(batch@runs[[j]]@plates[[k]]@mdata.unknown)&&dim(batch@runs[[j]]@plates[[k]]@mdata.unknown)[1]!=0)
				{
					dfm<-batch@runs[[j]]@plates[[k]]@mdata.unknown;
				}else {
					dfm<-aggregate(unknown[,"OD"], FUN=mean,by=list(unknown$ID))
					colnames(dfm)<-c("ID","OD");
				}
				colnames(dfm)[colnames(dfm)=="OD"]<-"OD_avg";
				for(q in 1:nRep)
				{
					dfm<-cbind(dfm, raw=NaN)
				}
				#extrCol.averageOD<-0;
				#if(is.null(batch@runs[[j]]@plates[[k]]@mdata.unknown)||dim(batch@runs[[j]]@plates[[k]]@mdata.unknown)[1]==0){
				#	dfm<-cbind(dfm, OD=NaN);
				#	extrCol.averageOD<-1;
				#}
				#rownames(dfm)<-dfm$ID;
				#dfm<-dfm[ids,]
				for(p in 1:length(ids))
				{
						smp<-unknown[unknown$ID==ids[p],]$OD
						dfm[dfm$ID==ids[p],c(dim(dfm)[2]-nRep+c(1:length(smp)))]<-smp;
						#if(is.null(batch@runs[[j]]@plates[[k]]@mdata.unknown)||dim(batch@runs[[j]]@plates[[k]]@mdata.unknown)[1]==0){
						#	dfm[dfm$ID==ids[p],"OD"]<-mean(smp);
						#}
				}
				#nRep<-nRep+extrCol.averageOD;
				dfm<-dfm[,c(1,dim(dfm)[2]-c(nRep:1)+1,c(2:(dim(dfm)[2]-nRep))) ]
			
				#now save the data
				suppressWarnings(write.table(dfm,file=file.conn,append=T,sep="\t", row.names = F))
				write("\n",file.conn, append=T);
			}
		}
	}
	close(file.conn);
}