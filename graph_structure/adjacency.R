source('header_global.R')
source('preprocess/header_preprocess.R')
library(plyr)

load(paste0(PATH_DATA, 'ABIDE_50002_matrix_2015-12-07.RData'))
load(paste0(PATH_DATA, 'template_2015-12-07.RData'))
MNI = load.nifti(paste0(PATH_DATA, 'MNI152_T1_2mm_brain_symmetric.nii.gz'))

dimen = dim(MNI)
mask = template$mask

loc = sapply(mask, .convert.2Dto3Dloc, dimen = dimen, simplify = 'matrix')

link.cross1 = function(ijk, dimen)
{
  i = ijk[1]; j = ijk[2]; k = ijk[3]
  
  a = matrix(c(i  , j  , k  ,
               i  , j  , k  ,
               i  , j  , k  ,
               i+1, j  , k  ,
               i  , j+1, k  ,
               i  , j  , k+1), nrow = 3)
  apply(a, 2, .convert.3Dto2Dloc, dimen = dimen)
}
link.cross2 = function(ijk, dimen)
{
  i = ijk[1]; j = ijk[2]; k = ijk[3]

  b = matrix(c(i  , j+1, k+1,
               i+1, j  , k+1,
               i+1, j+1, k  ,
               i  , j+1, k  ,
               i  , j  , k+1,
               i+1, j  , k  ), nrow = 3)
  apply(b, 2, .convert.3Dto2Dloc, dimen = dimen)
}

link1 = as.vector(apply(loc, 2, link.cross1, dimen = dimen))
link2 = as.vector(apply(loc, 2, link.cross2, dimen = dimen))
valid.link = link1 %in% mask & link2 %in% mask
link = cbind(link1, link2)[valid.link,]
rm(link1, link2, valid.link)

link = plyr::mapvalues(link, mask, 1:length(mask))

source('energy_parcellation/header_energy.R')
weights = compute.edgeWeights(dat$dat, NULL, dcor, link, verbose = T, save = F)
