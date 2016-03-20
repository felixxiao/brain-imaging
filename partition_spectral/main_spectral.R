source('partition_spectral/header_spectral.R')
source('criteria/header_criteria.R')

load(paste0(PATH_DATA, 'results/edges_ABIDE50002_nz.RData'))
part_nz.rc = partition.spectral.recursive(edges_nz, max.depth = 7)
part_nz.mt = partition.spectral.multiway(edges_nz, 128)

save(part_nz.rc, part_nz.mt,
     file = paste0(PATH_DATA, 'results/part_spec_128_ABIDE50002_nz.RData'))
