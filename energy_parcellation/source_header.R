suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(energy))
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(assertthat))

set.seed(10)
DATE = Sys.Date()

#path will vary from user to user
#where are you reading the data from?
PATH_DATA = "/home/smile/klsix/felix_senior_thesis_2015-16/data/"
#PATH_DATA = "C:/Users/Felix/Dropbox/Felix_Kevin_Han-seniorthesis2015-16/data/"

#where are the saving the outputs to?
PATH_SAVE = "/home/smile/klsix/felix_senior_thesis_2015-16/results/"
#PATH_SAVE = PATH_DATA

#test PATHs
#if PATHs don't end with a backslash, fail and purposely don't load the
# scripts.
for(i in c(PATH_DATA, PATH_SAVE)){
  tmp = rev(unlist(strsplit(i, "")))[1]
  assert_that(tmp == "/")
}


source("source_graph.R")
