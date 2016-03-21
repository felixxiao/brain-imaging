source('partition_edgecontract/header_contraction.R')

brain = 'ABIDE50002_nz'
load(paste0(PATH_DATA, 'results/edges_', brain, '.RData'))

start.time = Sys.time()
cg = partition.contractedge(5000, edges, T, brain)
print(Sys.time() - start.time)
