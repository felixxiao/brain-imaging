#going to be writing code pretty messily

setwd("C:/Users/UikosPC/Dropbox/Collaboration and People/Felix_Kevin/Felix_Kevin_Han-seniorthesis2015-16/data")
library(oro.nifti)
library(rgl)
library(misc3d)

template <- readNIfTI("MNI152_T1_2mm_brain.nii.gz")
template <- template@.Data
dtemp <- dim(template)


target.val = c(9032, 9140) #cerebellum_4_5_R and Vermis_7

aal = readNIfTI("AAL_MNI_2mm.nii")
aal = aal@.Data

par3d("windowRect" = c(0,23,1366,728))
contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3],
          level = 5000, alpha = .2, 
          #screen = list(x = 90, y = 90, z = 90), 
          light = c(0,0,-1))

aal1 = aal
aal1[which(aal != target.val[1])] = 0
contour3d(aal1, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], add = TRUE,
          level = 9032, color = rgb(219, 51, 64, max = 255))

aal2 = aal
aal2[which(aal != target.val[2])] = 0
contour3d(aal2, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], add = TRUE,
          level = 9032, color = rgb(45, 151, 227, max = 255))

setwd("C:/Users/UikosPC/Dropbox/Random/DUMP_figures")
rgl.snapshot("aal_backview.png")


par3d("windowRect" = c(0,23,1366,728))
par3d("userMatrix" = rotationMatrix(0, 1, 0, 0))
contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3],
          level = 5000, alpha = .2, 
          screen = list(x = 180, y = 0, z = 0), 
          light = c(0,0,-1))

aal1 = aal
aal1[which(aal != target.val[1])] = 0
contour3d(aal1, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], add = TRUE,
          level = 9032, color = rgb(219, 51, 64, max = 255))

aal2 = aal
aal2[which(aal != target.val[2])] = 0
contour3d(aal2, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], add = TRUE,
          level = 9032, color = rgb(45, 151, 227, max = 255))

setwd("C:/Users/UikosPC/Dropbox/Random/DUMP_figures")
rgl.snapshot("aal_topview.png")


par3d("windowRect" = c(0,23,1366,728))
rot.mat = matrix(0,4,4)
rot.mat[3,1] = rot.mat[1,2] = rot.mat[2,3] = rot.mat[4,4] = 1
par3d("userMatrix" = rot.mat)
contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3],
          level = 5000, alpha = .2, 
          screen = list(x = 180, y = 0, z = 0), 
          light = c(0,0,-1))

aal1 = aal
aal1[which(aal != target.val[1])] = 0
contour3d(aal1, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], add = TRUE,
          level = 9032, color = rgb(219, 51, 64, max = 255),light = c(-1,0,0))

aal2 = aal
aal2[which(aal != target.val[2])] = 0
contour3d(aal2, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], add = TRUE,
          level = 9032, color = rgb(45, 151, 227, max = 255))

setwd("C:/Users/UikosPC/Dropbox/Random/DUMP_figures")
rgl.snapshot("aal_topview.png")





#######
#goofing around:
contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3],
          level = 5000, alpha = .2, add=TRUE)


contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3],
          level = c(5000,8000), color = c("white", "red"), alpha = .2, add=TRUE)

template2 = template
template2[which(template2 != 0)] = 100

contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3],
          level = 100, alpha = .2, add=TRUE)