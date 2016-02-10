source('energy_parcellation/header_energy.R')

load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData')) #loads dat
#load(paste0(PATH_DATA, 'template_2015-12-07.RData')) #loads template
#load(paste0(PATH_DATA, 'edges_2016-01-14.RData')) #loads edges

rm.zero = preprocess.remove_zero_vertices(dat$dat)
dat = rm.zero$data
map = rm.zero$map
rm(rm.zero)

#edge.mat = preprocess.map_vertices(edges$edge.mat, map, inverse = T)
#edges = compute.edgeWeights(dat, NULL, edge.mat = edge.mat)

load(paste0(PATH_DATA, 'edges_2016-01-29.RData'))

part100 = partition.addedge.unconstrained(edges, 100)
table(part100)
table(table(part100))

