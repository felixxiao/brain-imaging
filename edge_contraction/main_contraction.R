source('edge_contraction/header_contraction.R')

load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

start.time = Sys.time()
partition.contractedge(edges, 230000)
print(Sys.time() - start.time)
