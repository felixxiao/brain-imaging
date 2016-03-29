source('package.R')
suppressMessages(.ipak(c("testthat", "assertthat", 'oro.nifti')))

USER = paste(Sys.getenv(c('USER', 'USERNAME')), collapse = '')

if (USER == "kevin" | USER == "klsix") {
  PATH_DATA = "/home/smile/klsix/felix_senior_thesis_2015-16/data/"
  PATH_SAVE = "/home/smile/klsix/felix_senior_thesis_2015-16/results/"
} else if (USER == "Felix") {
  PATH_DATA = "C:/Users/Felix/Dropbox/Felix_Kevin_Han-seniorthesis2015-16/data/"
  PATH_SAVE = 'C:/Users/Felix/Dropbox/Felix_Kevin_Han-seniorthesis2015-16/figs/'
} else if (USER == 'feixiao') {
  PATH_DATA = '/home/smile/feixiao/Thesis_data/'
  PATH_SAVE = PATH_DATA
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
  frame_files = lapply(sys.frames(), function(x) x$ofile)
  frame_files = Filter(Negate(is.null), frame_files)
  script.dir = dirname(frame_files[[length(frame_files)]])

  script.dir = normalizePath(script.dir)

  lapply(grep(glob2rx('^source_*.R'), dir(script.dir), value = T), function(x){
    source(paste0(script.dir, '/', x))
  })

  invisible()
}

load.library = function(package.names){
  suppressMessages(.ipak(package.names))

  invisible()
}

load.nifti = function(path)
{
  dat = readNIfTI(path)
  dat = dat@.Data

  dat
}
