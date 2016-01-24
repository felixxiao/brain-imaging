source('edge_contraction/header_contraction.R')

load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

start.time = Sys.time()
part = partition.contractedge(edges, 10)
print(Sys.time() - start.time)

save(part, file = paste0(PATH_SAVE, 'partition_contract_10_', DATE, '.RData'))
