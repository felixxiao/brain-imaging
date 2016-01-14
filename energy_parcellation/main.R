rm(list=ls())

load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
load(paste0(PATH_DATA, 'template_2015-12-07.RData'))
source("source_header.R")

construct.graph(dat$dat, template$neighbor.list, 20)
