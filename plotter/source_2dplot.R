#plot slices of the 3d cube
plot.partition2D <- function(partition, templateMRI, mask, filename,
  num.slices = 12, view = "sagittal") {
  assert_that(length(dim(templateMRI)) == 3)
  assert_that(length(partition) == length(mask))
  assert_that(length(partition) <= prod(dim(templateMRI)))
  assert_that(max(mask) <= prod(dim(templateMRI)))  

  #convert the partition into consecutive numerics
  vec = numeric(length(partition))
  lev = levels(partition)

  for(i in 1:length(lev)) { vec[which(partition == lev[i])] = i }

  #determine the slices to be plotted
  res = .compute.slices(templateMRI, view, num.slices)

  #replace the values in the templateMRI with those in partition
  templateMRI[mask] = vec
  image.slices = .extract.slices(templateMRI, res$slice.idx, res$dim.idx)

  #now plot the partition
  col.vec = c("black", sample(rainbow(length(lev))))
  
  res = .compute.plotLayout(num.slices)

  png(file = paste0(PATH_SAVE, "partition_", filename, "_" ,
   view, "_", DATE, ".png"), 
   height = 4, width = 6, units = "in", res = 600)

  par(mfrow = c(res$num.row, res$num.col), mar = rep(0.2, 4), bg = "black")

  for(i in 1:length(image.slices)){
    image(image.slices[[i]], col = col.vec,
     breaks = seq(-0.5, 20.5, length.out = length(col.vec)+1), 
     zlim = c(0, max(lev)), bty = "n", xaxt = "n", yaxt = "n",
     asp = ncol(image.slices[[1]])/nrow(image.slices[[1]]))
  }

  dev.off()
  invisible()
} 

.compute.plotLayout <- function(num.plots) {
  num.row = ceiling(sqrt(num.plots/2))
  num.col = ceiling(num.plots / num.row)

  list(num.row = num.row, num.col = num.col)
}


#given a 3D image (img), split img along the dim.idx
#  direction, and output only the slice.indices
.extract.slices <- function(img, slice.idx, dim.idx){
  img.slices = alply(img, dim.idx)
  img.slices[slice.idx]
}


#given an image, find out the dimensional to slice along and the
# slices along that dimension
.compute.slices <- function(img, view, num.slices = 12, offset = 5) {  
  assert_that(view %in% VIEWS & length(view) == 1)
  
  dim.idx = which(VIEWS %in% view)  
  num.slices = min(num.slices, dim(img)[dim.idx])

  if(offset > dim(img)[dim.idx]) offset = 0

  #determine which dimension indicies have actual "image data" in them
  tmp = apply(img, dim.idx, function(x){sum(abs(x))})
  sort.tmp = sort(which(tmp>0), decreasing = FALSE)

  min.slice = sort.tmp[1] + offset
  max.slice = rev(sort.tmp)[1] - offset
  
  slice.idx = round(seq(min.slice, max.slice, length.out = num.slices))
  
  list(slice.idx = slice.idx, dim.idx = dim.idx)
}
