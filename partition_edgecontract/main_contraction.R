source('partition_edgecontract/header_contraction.R')
source('criteria/header_criteria.R')

brain = 'ABIDE50002_nz'
load(paste0(PATH_DATA, 'results/edges_', brain, '.RData'))
load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
dat = dat$dat[,map_nz]
edges_nz = preprocess.validate.edges(edges_nz)

num_components = c(300, 116)
alpha = c(3, 6, 10, 6, 6)
beta  = c(1, 1,  1, 2, 4)

for (i in 1:length(alpha))
{
  partitions = partition.contractedge.general(num_components, 
                                              alpha[i], beta[i], edges_nz, NULL,
                                              length(map_nz))$partitions
  for (k in 1:length(num_components))
  {
    cat(num_components[k], ' components\n')
    cat('Adjac: ', criterion.adjacent_pairwise_ecor(edges_nz, partitions[[k]])$total.mean, '\n')
    #  cat('Multi: ', criterion.multi_boundary_ecor(dat, partitions[[k]]), '\n')
    cat('Jaggd: ', criterion.jaggedness(edges_nz, partitions[[k]]), '\n')
    cat('Balan: ', criterion.balance(partitions[[k]]), '\n')
    cat('\n')
  }
}
