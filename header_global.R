suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(testthat))

#USER = "kevin"
USER = "felix"

if (USER == "kevin") {
  PATH_DATA = "/home/smile/klsix/felix_senior_thesis_2015-16/data/"
  PATH_SAVE = "/home/smile/klsix/felix_senior_thesis_2015-16/results/"
} else if (USER == "felix") {
  PATH_DATA = "C:/Users/Felix/Dropbox/Felix_Kevin_Han-seniorthesis2015-16/data/"
  PATH_SAVE = 'C:/Users/Felix/Dropbox/Felix_Kevin_Han-seniorthesis2015-16/figs/'
}

# verify PATHs end in '/'
for (p in c(PATH_DATA, PATH_SAVE))
	assert_that(substr(p, nchar(p), nchar(p)) == '/')
rm(p)

set.seed(10)
DATE  = Sys.Date()
VIEWS = c("sagittal", "coronal", "axial")

######################################################################

# assortment of functions for easy and universal loading of functions
# to be used by each folder

# load all the source files in current directory
load.source = function()
{
  lapply(grep("^source*", dir(), value = TRUE), source)
}

load.nifti = function(path)
{
  dat = readNIfTI(path)
  dat = dat@.Data
}
