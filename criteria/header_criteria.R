source('header_global.R')

# load all validation_criteria/source_* files
lapply(grep('^source_', dir('criteria'), value = T), source)
