source('preprocess/header_preprocess.R')
source('graph_structure/header_structure.R')
source('energy_parcellation/header_energy.R')

load(paste0(PATH_DATA, 'template_2015-12-07.RData'))
MNI = load.nifti(paste0(PATH_DATA, 'MNI152_T1_2mm_brain_symmetric.nii.gz'))

dimen = dim(MNI)
mask = template$mask
rm(MNI, template)

loc = sapply(mask, .convert.2Dto3Dloc, dimen = dimen, simplify = 'matrix')

edge.mat = list()
for (dist_sqrd in 1:16)
{
  time = Sys.time()
  edge.mat[[dist_sqrd]] = edgemat(dist_sqrd, loc, dimen, mask)
  print(paste(dist_sqrd, ':', Sys.time() - time))
  save(edge.mat, file = paste0(PATH_SAVE, 'edgemat_dist_1-12.RData'))
}

#load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
#weights = compute.edgeWeights(dat$dat, NULL, dcor, edge.mat2, verbose = T, save = F)
#weights = weights$energy.vec
