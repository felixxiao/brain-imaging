rm(list=ls())

source("source_header.R")
load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData')) #loads dat
load(paste0(PATH_DATA, 'template_2015-12-07.RData')) #loads template
load(paste0(PATH_SAVE, 'edges_2016-01-14.RData')) #loads edges

g = construct.graph(NULL, NULL, edges, component.num = 20)
save(g, file = paste0(PATH_SAVE, "graph_", DATE, ".RData"))

comp = components(g)
unique(dat$dat[,which(comp$membership %in% 2:20)])
zeros = which(apply(dat$dat, 2, function(x) all(unique(x) == 0)))
