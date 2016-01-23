rm(list = ls())

source("header_plotter.R")

load(paste0(PATH_DATA, 'template_2015-12-07.RData')) 
MNI = load.nifti(paste0(PATH_DATA, 'MNI152_T1_2mm_brain_symmetric.nii.gz'))

#create some random 
partition = as.factor(sample(1:20, length(template$mask), replace = TRUE))
plot.partition2D(partition, MNI, template$mask, "test")
plot.partition2D(partition, MNI, template$mask, "test", view = VIEWS[2])
plot.partition2D(partition, MNI, template$mask, "test", view = VIEWS[3])



