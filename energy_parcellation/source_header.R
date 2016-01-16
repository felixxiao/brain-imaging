suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(energy))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(assertthat))



source("source_graph.R")

set.seed(10)
DATE = Sys.Date()

#path will vary from user to user
#where are you reading the data from?
#PATH_DATA = "~/fmri_script_test/20151206_felix_output/"
PATH_DATA = "C:/Users/Felix/Dropbox/Felix_Kevin_Han-seniorthesis2015-16/data/"

#where are the saving the outputs to?
#PATH_SAVE = "~/fmri_script_test/20151206_felix_output/"
PATH_SAVE = PATH_DATA

