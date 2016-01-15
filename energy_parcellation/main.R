rm(list=ls())

source("source_header.R")
load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
load(paste0(PATH_DATA, 'template_2015-12-07.RData'))

g = construct.graph(dat$dat, template$neighbor.list, 20)
save(g, file = paste0(PATH_SAVE, "graph_", DATE, ".RData"))
