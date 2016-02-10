# for all instances in dimension 4, compress along dimensions 1...3
# by factor, taking the mean of the factor^3 cells in each block
compress.Array4D = function(array.4d, factor)
{
  .validate.Array4D(array.4d)
  assert_that(factor > 1)

  A = array.4d
  
  dim.comp = dim(A)
  dim.comp[1:3] = dim.comp[1:3] %/% factor

  idx = seq(1, by = factor, length.out = dim.comp[1])
  B = A[idx,,,]
  for (i in 2:factor)
    B = B + A[idx + i - 1,,,]

  idx = seq(1, by = factor, length.out = dim.comp[2])
  A = B[,idx,,]
  for (i in 2:factor)
    A = A + B[,idx + i - 1,,]

  idx = seq(1, by = factor, length.out = dim.comp[3])
  B = A[,,idx,]
  for (i in 2:factor)
    B = B + A[,,idx + i - 1,]

  B / factor^3
}

compress.mask = function(mask, dimen, factor)
{
  .validate.mask(mask, dimen)
  mask.3d = lapply(mask, .convert.2Dto3Dloc, dimen = dimen)

  dimen = dimen %/% factor
  mask.comp = numeric(0)
  for (i in 1:length(mask.3d))
    if (all(mask.3d[[i]] %% factor == 0))
      mask.comp = c(mask.comp,
                    .convert.3Dto2Dloc(mask.3d[[i]] / factor, dimen))
  mask.comp
}

convert.Array4D.Matrix = function(array.4d, mask)
{
  .validate.mask(mask, dim(array.4d)[1:3])
  .validate.Array4D(array.4d)

  t(apply(array.4d, 4, function(x) x[mask]))
}

convert.mask.template = function(mask, dimen)
{
  .validate.mask(mask, dimen)
  template = array(0, dim = dimen)
  template[mask] = 1
  template
}

compress.data = function(array.4d, mask, factor = 2, verbose = T)
{
  dimen = dim(array.4d)[1:3]
  .validate.mask(mask, dimen)
  .validate.Array4D(array.4d)

  if (verbose) cat('Compressing 4-D array\n')
  array.4d = compress.Array4D(array.4d, factor)
  
  if (verbose) cat('Creating new mask\n')
  mask = compress.mask(mask, dimen, factor)

  if (verbose) cat('Converting 4-D array to matrix\n')
  mat = convert.Array4D.Matrix(array.4d, mask)
  
  dimen = dim(array.4d)[1:3]

  if (verbose) cat('Creating new template\n')
  template = convert.mask.template(mask, dimen)

  if (verbose) cat('Extracting neighbors\n')
  neighbors = extract.neighbors(template, verbose = F)$neighbor.list

  list(mat = mat, mask = mask, template = template,
       neighbors = neighbors)
}

.validate.mask = function(mask, dimen)
{
  assert_that(is.vector(mask) & is.numeric(mask))
  assert_that(max(mask) <= prod(dimen))
}

.validate.Array4D = function(array.4d)
{
  assert_that(class(array.4d) == 'array')
  assert_that(length(dim(array.4d)) == 4)
}

# given the MNI standard template, find out all the neighbors
# return two things: the mask and list of neighbors
# the pattern dictates how neighbors are defined.
extract.neighbors = function(template, pattern = .cross.enumerate(), 
                             tmpsave = NA, verbose = T)
{
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

# convert a location (single matrix index) into 3D coordinates
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

# convert a location (3D by coordinates) into an index for a matrix
.convert.3Dto2Dloc = function(loc, dimen)
{
  assert_that(length(loc)==3 & is.numeric(loc))
  assert_that(length(dimen)==3 & is.numeric(dimen))
  assert_that(all(loc<=dimen))
  assert_that(all(loc>0))

  loc[1]+dimen[1]*(loc[2]-1)+(dimen[1]*dimen[2])*(loc[3]-1)
}
