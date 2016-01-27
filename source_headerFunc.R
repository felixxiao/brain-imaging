suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(testthat))

#set global variables
#path will vary from user to user
#where are you reading the data from?
USER = "kevin"
#USER = "felix"

if(USER == "kevin") {
  PATH_DATA = "/home/smile/klsix/felix_senior_thesis_2015-16/data/"
  PATH_SAVE = "/home/smile/klsix/felix_senior_thesis_2015-16/results/"
} else if(USER == "felix") {
  PATH_DATA = "C:/Users/Felix/Dropbox/Felix_Kevin_Han-seniorthesis2015-16/data/"
  PATH_SAVE = PATH_DATA
}

#test PATHs
#if PATHs don't end with a backslash, fail and purposely don't load the
# scripts.
for(i in c(PATH_DATA, PATH_SAVE)){
  tmp = rev(unlist(strsplit(i, "")))[1]
  assert_that(tmp == "/")
}

rm(list = c("i", "tmp"))

set.seed(10)
DATE  = Sys.Date()
VIEWS = c("sagittal", "coronal", "axial")

##################

#assortment of functions for easy and universal loading of functions
#to be used by each folder

#load all the source files in the directory
load.source <- function() {
  lapply(grep("^source*", dir(), value = TRUE), source)
}

load.nifti <- function(path){
  dat = readNIfTI(path)
  dat = dat@.Data
}



