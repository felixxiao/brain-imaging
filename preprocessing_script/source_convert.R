#given fmri data (3D) and a matrix mask, reformat it into a matrix
extract.data <- function(dat, mask, verbose = TRUE){
  assert_that(is.numeric(dat) & length(dim(dat))==4)
  assert_that(is.vector(mask) & is.numeric(mask))
  assert_that(length(mask) <= length(dat))

  dimen = dim(dat)
  assert_that(max(mask) <= prod(dimen[1:3]))

  mat = matrix(NA, dimen[4], length(mask))
  for(j in 1:dimen[4]){
    mat[j,] = as.numeric(dat[,,,j])[mask]

    if(j%%floor(dimen[4]/10)==0 & verbose) cat('*')
  }

  mat
}

#given the MNI standard template, find out all the neighbors
#return two things: the mask and list of neighbors
#the pattern dictates how neighbors are defined.
extract.neighbors <- function(template, pattern = .cross.enumerate(), 
 tmpsave = NA, verbose = TRUE){

  assert_that(is.numeric(template) & length(dim(template))==3)
  assert_that(is.matrix(pattern) & ncol(pattern)==3)

  dimen = dim(template)
  mask = which(template != 0)
  neighbor.list = list(mask)

  for(j in 1:length(mask)){
    #convert mask index into pixel-location
    loc = .convert.2Dto3Dloc(mask[j], dimen)
		
    #apply pattern
    neigh = t(apply(pattern,1,function(s, loc){s+loc}, loc=loc))

    #remove the invalid neighbors
    rmv = unique(c(which(neigh<=0,arr.ind=TRUE)[,1],
     which(neigh[,1]>dimen[1]), which(neigh[,2]>dimen[2]),
     which(neigh[,3]>dimen[3])))
    if(length(rmv)>0) neigh = neigh[-rmv,]

    #convert pattern back to idx and see if it's in mask
    idx = apply(neigh, 1, .convert.3Dto2Dloc, dimen=dimen)
    #see if the index number is in mask
    idx = idx[which(idx %in% mask)]
    #convert idx into a column number
    colidx = which(mask %in% idx)

    #make sure there were neighbors found
    assert_that(length(colidx) == length(idx))
    if(length(colidx)==0) {
      stop(paste0("ERROR AT INDEX ",j," WHERE NO NEIGHBORS FOUND!"))
    }
		
    neighbor.list[[j]] = colidx
    assert_that(max(colidx) <= length(mask))
    assert_that(min(colidx) > 0)

    if(length(mask)>1000 & j %% floor(length(mask)/1000)==0) {
      if(verbose) cat('*')
      if(!is.na(tmpsave)) save(neighbor.list, file=tmpsave)
    }
  }

  assert_that(length(neighbor.list)==length(mask))
  list(neighbor.list = neighbor.list, mask = mask)
}


.cube.enumerate <- function(){
  vec = c( 0, 0, 1,   0, 1, 0,  1, 0, 0,
           0, 0,-1,   0,-1, 0, -1, 0, 0,
           0, 1, 1,   1, 0, 1,  1, 1, 0,
           0,-1, 1,  -1, 0, 1, -1, 1, 0,
           0, 1,-1,   1, 0,-1,  1,-1, 0,
           0,-1,-1,  -1, 0,-1, -1,-1, 0,
           1, 1, 1,  -1,-1,-1,
           1,-1,-1,  -1, 1,-1, -1,-1, 1,
           1, 1,-1,   1,-1, 1, -1, 1, 1)

  mat = matrix(vec,ncol=3,byrow=TRUE)
  mat
}

.cross.enumerate <- function(){
  vec = c( 0, 0, 1,   0, 1, 0,  1, 0, 0,
           0, 0,-1,   0,-1, 0, -1, 0, 0)

  mat = matrix(vec,ncol=3,byrow=TRUE)
  mat
}

#convert a location (in matrix represented by an index) 
#  into 3D coordinates
.convert.2Dto3Dloc <- function(idx, dimen){
  assert_that(is.numeric(idx))
  assert_that(length(dimen)==3 & is.numeric(dimen))
  assert_that(idx <= prod(dimen))

  z = ceiling(idx / (dimen[1]*dimen[2]))

  tmp = idx %% (dimen[1]*dimen[2])
  if(tmp==0) tmp = dimen[1]*dimen[2]
  y = ceiling(tmp / dimen[1])

  x = tmp %% dimen[1]
  if(x==0) x = dimen[1]

  c(x,y,z)
}

#convert a location (3D by coordinates) into an index for
#  a matrix
.convert.3Dto2Dloc <- function(loc, dimen){
  assert_that(length(loc)==3 & is.numeric(loc))
  assert_that(length(dimen)==3 & is.numeric(dimen))
  assert_that(all(loc<=dimen))
  assert_that(all(loc>0))

  loc[1]+dimen[1]*(loc[2]-1)+(dimen[1]*dimen[2])*(loc[3]-1)
}
