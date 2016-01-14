suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(energy))

source("source_graph.R")

set.seed(10)
DATE = Sys.Date()

#path will vary from user to user
PATH_DATA = "~/fmri_script_test/20151206_felix_output/"
#PATH_DATA = "C:/Users/Felix/Dropbox/Felix_Kevin_Han-seniorthesis2015-16/data"/

PATH_SAVE = "~/DUMP/"
