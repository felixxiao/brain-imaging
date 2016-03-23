source('partition_edgecontract/header_contraction.R')
source('criteria/header_criteria.R')

brain = 'ABIDE50002_nz'
load(paste0(PATH_DATA, 'results/edges_', brain, '.RData'))
load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
dat = dat$dat[,map_nz]
edges_nz = preprocess.validate.edges(edges_nz)

num_components = c(500, 400, 300, 250, 200, 150, 116, 100)
partitions = partition.contractedge.python(num_components, edges_nz)$partitions

for (k in 1:length(num_components))
{
  cat(num_components[k], ' components\n')
  cat('Adjac: ', criterion.adjacent_pairwise_ecor(edges_nz, partitions[[k]])$total.mean, '\n')
#  cat('Multi: ', criterion.multi_boundary_ecor(dat, partitions[[k]]), '\n')
  cat('Jaggd: ', criterion.jaggedness(edges_nz, partitions[[k]]), '\n')
  cat('Balan: ', criterion.balance(partitions[[k]]), '\n')
  cat('\n')
}
