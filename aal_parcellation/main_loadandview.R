filename = '~/felix_senior_thesis_2015-16/data/AAL_MNI_2mm.nii'
template = '~/felix_senior_thesis_2015-16/data/template_2015-12-07.RData'
mni.file = '~/felix_senior_thesis_2015-16/data/MNI152_T1_2mm_brain_symmetric.nii.gz'

load(template)
adj.list = template$neighbor.list
mask = template$mask

source("energy_parcellation/header_energy.R")
source("aal_parcellation/source_read.R")
library(assertthat)
library(oro.nifti)

aal.partition = read.aal(filename, adj.list, mask)
assert_that(length(aal.partition) == length(mask))

library(plyr)
source("plotter/source_2dplot.R")

mni = readNIfTI(mni.file)
mni = mni@.Data

set.seed(10)
DATE  = Sys.Date()
VIEWS = c("sagittal", "coronal", "axial")
PATH_SAVE = '~/felix_senior_thesis_2015-16/results/'

plot.partition2D(aal.partition, mni, mask, filename = 'aal',
  num.slices = 24, height = 8, width = 12)
plot.partition2D(aal.partition, mni, mask, filename = 'aal',
  num.slices = 24, view = VIEWS[3], height = 8, width = 12)

