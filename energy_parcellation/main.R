rm(list=ls())

source("source_header.R")
load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
load(paste0(PATH_DATA, 'template_2015-12-07.RData'))
load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

g = construct.graph(NULL, NULL, edges, component.num = 20)
save(g, file = paste0(PATH_SAVE, "graph_", DATE, ".RData"))

components = components(g)
unique(dat$dat[,which(components$membership %in% 2:20)])
zeros = which(apply(dat$dat, 2, function(x) all(unique(x) == 0)))
