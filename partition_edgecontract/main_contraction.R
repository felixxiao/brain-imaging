source('partition_edgecontract/header_contraction.R')
source('criteria/header_criteria.R')

brain = 'ABIDE50002_nz'
load(paste0(PATH_DATA, 'results/edges_', brain, '.RData'))

write.table(edges_nz$edge.mat, 'partition_edgecontract/edge_mat.csv',
            sep = ',', eol = '\n', row.names = F, col.names = F)
write.table(edges_nz$energy.vec, 'partition_edgecontract/weights.csv',
            eol = '\n', row.names = F, col.names = F)

# call Python

partitions = lapply(strsplit(readLines('partition_edgecontract/partition.csv'), ','), as.factor)
num_components = partitions[[1]]
partitions = partitions[2:length(partitions)]

load(paste0(PATH_DATA, 'results/edges_ABIDE50002_nz.RData'))
load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
dat = dat$dat[,map_nz]
dat = dat[, 1:ncol(dat) != 30351]

for (k in 2:length(num_components))
{
  cat(num_components[k], ' components\n')
  cat('Adjac: ', criterion.adjacent_pairwise_ecor(edges_nz, partitions[[k]])$total.mean, '\n')
  cat('Multi: ', criterion.multi_boundary_ecor(dat, partitions[[k]]), '\n')
  cat('S-V-r: ', criterion.surface_volume_ratio(edges_nz, partitions[[k]]), '\n')
  cat('Balan: ', criterion.balance(partitions[[k]]), '\n')
  cat('\n')
}
