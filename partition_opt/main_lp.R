############ Compress ABIDE 50002 by 4x4x4 ##########################
source('preprocess/header_preprocess.R')
load(paste0(PATH_DATA, 'template_2015-12-07.RData'))
array.4d = load.nifti(paste0(PATH_DATA, 'func2mni.nii.gz'))

comp.4 = compress.data(array.4d, template$mask, factor = 4)
rm(array.4d, template)

source('energy_parcellation/header_energy.R')
comp.4$edges = compute.edgeWeights(comp.4$mat, comp.4$neighbors)

save(comp.4, file = paste0(PATH_DATA, 'ABIDE_50002_compress_4.RData'))
#####################################################################
load(paste0(PATH_DATA, 'ABIDE_50002_compress_4.RData'))

source('partition_opt/header_opt.R')

L = compute.laplacian(comp.4$edges)

write.ampl_data(L, 3, 'data.ampl')

# call main.ampl
# compress_4 n = 3666

#parcel = compute.factor_01.julia('partition_opt/Z.csv', 20)
X = read.table('partition_opt/X.csv', sep = ',')
parcel = rep(NA, times = nrow(X))
for (j in 1:ncol(X)) parcel[X[,j] == 1] = j
parcel = as.factor(parcel)

#####################################################################
source('plotter/header_plotter.R')

plot.partition2D(parcel, comp.4$template, comp.4$mask, filename = 'comp_4_lp')
