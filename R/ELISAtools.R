###ELISA tools, an R software package for analyzing ELISA data and 
#	doing batch effect correction.
#	By Feng @ BU 07/01/2018 Boston University (ffeng@bu.edu)
#
#  Version 0.1.5
#  For now, 1)it will only fit either the five- or four-parameter logistic model
#  2)do analysis to predicate based on the calibration the unknown concentration
#  3)if there are batches, will fit to estimate the shift factor and do correction/normalization 
#
#### Developed By Feng @ BU 2018. All right reserved###

#Below are roxygen2 comments

##########Dependency###########
##  library(MASS)
####################
#this is where the files to be include if necessary
#' @include ELISAtools_IO.R
#'
#'
#' @title ELISA data analysis with batch correction
#'
#' @description An R package to run ELISA data analysis 
#'		with the ability to do batch correction/normalization 
#'		
#'
#' @details This package is developed to run analysis of ELISA data. 
#'		First, the calibration data are used to fit either the five- or
#'		four-parameter logistic model. Then the fitted model is
#'		used to predict the concentrations of unknown samples.
#'		If the batches of calibration data exist, 
#'		the correction/normalization could be done. The corrected
#'		calibration curve are then used for predication.
#'
#' 		Please refer to the vignettes to see details.
#' @references Feng, et al 2018 \doi{10.1101/483800}		
"_PACKAGE"