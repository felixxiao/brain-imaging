source('preprocess/header_preprocess.R')
source('graph_structure/header_structure.R')
source('energy_parcellation/header_energy.R')

load(paste0(PATH_DATA, 'template_2015-12-07.RData'))
MNI = load.nifti(paste0(PATH_DATA, 'MNI152_T1_2mm_brain_symmetric.nii.gz'))

dimen = dim(MNI)
mask = template$mask
rm(MNI, template)

loc = sapply(mask, .convert.2Dto3Dloc, dimen = dimen, simplify = 'matrix')

time = Sys.time()
edge.mat2 = edgemat(2, loc, dimen, mask)
Sys.time() - time

load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
weights = compute.edgeWeights(dat$dat, NULL, dcor, edge.mat2, verbose = T, save = F)
weights = weights$energy.vec

mean(weights)
hist(weights)

save(link, weights, file = paste0('graph_structure/2-3.RData'))
