setwd('~/GitHub/union_find')
source('header_unionfind.R')

load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

part1 = partition.addedge.uf(edges, max.size = 15000, min.size = 1000, PATH_SAVE)
#load(paste0(PATH_DATA, 'edge_connect_1000_15000_2016-01-19.RData'))

part2 = partition.addedge.uf(edges, max.size = 10000, min.size = 1000, PATH_SAVE)
