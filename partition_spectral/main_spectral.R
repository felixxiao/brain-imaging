source('partition_spectral/header_spectral.R')
#source('energy_parcellation/header_energy.R')

load(paste0(PATH_DATA, 'results/edges_ABIDE50002_nz.RData'))

part_nz.rc = partition.spectral.recursive(edges_nz, max.depth = 7)
part_nz.mt = partition.spectral.multiway(edges_nz, 128)

save(part_nz.rc, part_nz.mt,
     file = paste0(PATH_DATA, 'results/part_spec_128_ABIDE50002_nz.RData'))

part.rc = preprocess.map_partition(part.rc, map_nonzero, N)
part.mt = preprocess.map_partition(part.mt, map_nonzero, N)

source('criteria/header_criteria.R')
load(paste0(PATH_DATA, 'edges_2016-01-14.RData'))

criteria = c('Adjacent', 'Boundary')
methods  = c('Recursive', 'Multiway')
crit = matrix(NA, nrow = length(methods), ncol = length(criteria),
              dimnames=  list(methods, criteria))

crit['Recursive', 'Adjacent'] = criterion.adjacent_pairwise_ecor(edges, part.rc)
bou = criterion.boundary_pairwise_ecor(edges, part.rc)

source('plotter/header_plotter.R')
MNI = load.nifti(paste0(PATH_DATA, 'MNI152_T1_2mm_brain_symmetric.nii.gz'))
