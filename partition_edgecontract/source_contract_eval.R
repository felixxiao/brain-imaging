source('partition_edgecontract/header_contraction.R')
source('criteria/header_criteria.R')

load(paste0(PATH_DATA, 'results/contractible_graph_1000_2016-02-14.RData'))
load(paste0(PATH_DATA, 'results/edges_ABIDE50002_nz.RData'))
load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
dat = dat$dat[,map_nz]
rm(edges_nz, map_nz)

K = c(250, 200, 150, 100)
for (k in K)
{
  cg = partition.contractedge(k, cg, T)
  part = cg$get_vertex_components()
  cat(k, ' components\n')
#  cat('Adjac: ', criterion.adjacent_pairwise_ecor(edges_nz, part)$total.mean, '\n')
  cat('Multi: ', criterion.multi_boundary_ecor(dat, part), '\n')
#  cat('S-V-r: ', criterion.surface_volume_ratio(edges_nz, part), '\n')
#  cat('Balan: ', criterion.balance(part), '\n')
#  cat('\n')
}

"
which(part == which(as.vector(table(part)) == 1))

# recover original brain parcellation
source('energy_parcellation/header_energy.R')
load(paste0(PATH_DATA, 'results/edges_ABIDE50002_nz.RData'))
load(paste0(PATH_DATA, 'results/edges_2016-01-14.RData'))
preprocess.assign_unmapped_vertices(map_nz, edges, N, part)

source('criteria/header_criteria.R')

criterion.balance(comp)
criterion.surface_volume_ratio(edges, comp)
"