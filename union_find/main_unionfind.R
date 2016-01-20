source('header_unionfind.R')

load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

part1 = partition.addedge.uf(edges, max.size = 15000, min.size = 1000, PATH_SAVE)
