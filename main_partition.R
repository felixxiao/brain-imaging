source('energy_parcellation/header_energy.R')

brain = 'ABIDE/50002'
k = 116

# get dat and edges for brain
load(paste0(brain, '-matrix.RData'))
load(paste0(brain, '-edges.RData'))
dat = dat$dat
N = ncol(dat)

# form dat_nz and edges_nz
out = preprocess.all(dat, edges)
dat_nz = out$dat_nz
edges_nz = out$edges_nz
map_nz = out$map_nz
rm(out)

# run partitioning algorithms
part = list()

source('partition_edgecontract/header_contraction.R')
part$ec.3.1 = partition.contractedge.general(k, 3, 1, edges_nz)$partitions[[1]]
part$ec.6.1 = partition.contractedge.general(k, 6, 1, edges_nz)$partitions[[1]]
part$ec.6.4 = partition.contractedge.general(k, 6, 4, edges_nz)$partitions[[1]]

source('partition_spectral/header_spectral.R')
part$sp = partition.spectral.multiway(edges_nz, k)
partition.contractedge.write_cg_json(edges_nz, part$sp, 'cg_sp.json')
part$sp.ec = partition.contractedge.general(k, 6, 4, cg.json = 'cg_sp.json')$partitions[[1]]

save(part, file = paste0(brain, '-partitions.RData'))

source('partition_nmf/header_nmf.R')
part$bf = partition.sym_bmf.python(edges_nz, k)
partition.contractedge.write_cg_json(edges_nz, part$bf, 'cg_bf.json')
part$bf.ec = partition.contractedge.general(k, 6, 4, cg.json = 'cg_bf.json')$partitions[[1]]

save(part, file = paste0(brain, '-partitions.RData'))
