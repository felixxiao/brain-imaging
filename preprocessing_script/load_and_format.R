rm(list=ls())
library(oro.nifti)
source("common_func.R")

subject.vec = c("50002")

standard = readNIfTI("~/fmri_script_test/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz")
dimen = dim(standard)

standard = as.numeric(standard@.Data)
mask = which(standard!= 0)

for(i in 1:length(subject.vec)){
	subj.name = paste("ABIDE_",subject.vec[i],"_resting_res2standard.nii.gz",sep="")
#	func = readNIfTI(subj.name)


#	print(paste("Start extraction for ",subject.vec[i],sep=""))
	#extract matrix
#	dat = func@.Data
#	dimen = dim(dat)
#	mat = matrix(NA,dimen[4],length(mask))
#	for(j in 1:dimen[4]){
#  		mat[j,] = as.numeric(dat[,,,j])[mask]
#
#  		if(j%%floor(dimen[4]/10)==0) cat('*')
#	}

#	save(mat, file=paste("processed_",subject.vec[i],"_matrix_FCP-Preprocess_07092015.RData",sep=""))

	print("Start enumerating neighbors")
	#list the neighboring indices
	neighbor.list = list(0)
	converter = 1:length(mask)
	pattern = sphere_enumerate(5)
	
	for(j in 1:length(mask)){
		#convert column index j into pixel-location
		loc = mask[j]
		z = ceiling(loc / (dimen[1]*dimen[2]))
    		tmp = loc %% (dimen[1]*dimen[2])
    		if(tmp==0) tmp = dimen[1]*dimen[2]
    		y = ceiling(tmp / dimen[1])
    		x = tmp %% dimen[1]
    		if(x==0) x = dimen[1]

		#apply pattern
		neigh = t(apply(pattern,1,function(s){s+c(x,y,z)}))

		#convert pattern back to idx and see if it's in mask
		tmp = unique(c(which(neigh<=0,arr.ind=TRUE)[,1],
                   which(neigh[,1]>dimen[1]), which(neigh[,2]>dimen[2]),
                   which(neigh[,3]>dimen[3])))
    		if(length(tmp)>0) neigh = neigh[-tmp,]
    		idx = apply(neigh,1,function(x){x[1]+dimen[1]*(x[2]-1)+(dimen[1]*dimen[2])*(x[3]-1)})
    		idx = idx[which(idx %in% mask)]
    		idx = converter[which(mask %in% idx)]

		if(length(idx)==0) {
			print(paste("ERROR AT INDEX ",j," WHERE NO NEIGHBORS FOUND",sep=""))
			stop("algorithm stopped!")
		}
		neighbor.list[[j]] = idx

		if(length(neighbor.list)!=j){
		 	print(paste("LENGTH ERROR.ERROR AT INDEX ",j," WHERE NO NEIGHBORS FOUND",sep=""))
                        stop("algorithm stopped!")
		}

		if(j %% floor(length(mask)/1000)==0) {
			cat('*')
			save(neighbor.list, file=paste("tmp_processed_",subject.vec[i],"_neighboring-list_07092015.RData",sep=""))
		}
	}

	save(neighbor.list, file=paste("processed_",subject.vec[i],"_neighboring-list_07082015.RData",sep=""))

}
