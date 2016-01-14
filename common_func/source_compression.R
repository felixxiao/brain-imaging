#compress the data according to func applied to every 2x2x2 block
compress.data <- function(dat, func = mean, verbose = TRUE){
  assert_that(is.numeric(dat))
  assert_that(length(dim(dat))==3 | length(dim(dat))==4)
  assert_that(length(dat)>0)

  dimen = dim(dat)
  halfdimen = dimen
  halfdimen[1:3] = ceiling(halfdimen[1:3]/2)
  compress.dat = array(NA, halfdimen)

  if(length(dimen)==3) dimen = c(dimen,1) #artifically add in a 4th dimen

  compressfunc <- function(i, time){
    #compute 3D location of the voxel in the compressed data
    c.pos = convert.2Dto3Dloc(i, halfdimen[1:3])
    
    #compute the 3D locations of the voxel in the full data
    lower = 2*c.pos-1
    upper = sapply(1:3, function(x){min(2*c.pos[x], dimen[x])})

    compress.dat[c.pos[1], c.pos[2], c.pos[3], time] <<- func(dat[lower[1]:upper[1], 
     lower[2]:upper[2], lower[3]:upper[3], time])

    invisible()
  }

  for(s in 1:dimen[4]){
    sapply(1:prod(halfdimen[1:3]), compressfunc, time=s)
    if(verbose & (dimen[4] < 10 || s %% floor(dimen[4]/10)==0)) cat('*')
  }

  compress.dat
}

