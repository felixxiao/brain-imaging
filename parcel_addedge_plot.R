source('plotter/header_plotter.R')

# this contains my partitions with size-constrained add-edge
load(paste0(PATH_DATA, 'parcel_addedge_eval_2016-01-23.RData'))

load(paste0(PATH_DATA, 'template_2015-12-07.RData'))
MNI = load.nifti(paste0(PATH_DATA, 'MNI152_T1_2mm_brain_symmetric.nii.gz'))

p = 2 # 7 partitions in all
plot.partition2D(parcel[[p]], MNI, template$mask, filename = p, view = VIEWS[1])
plot.partition2D(parcel[[p]], MNI, template$mask, filename = p, view = VIEWS[2])
plot.partition2D(parcel[[p]], MNI, template$mask, filename = p, view = VIEWS[3])

levels(parcel[[p]]) # number of components
table(parcel[[p]])  # size of each component
