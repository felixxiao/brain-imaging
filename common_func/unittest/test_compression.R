library(testthat)
setwd("..")

source("source_header.R")

test_that("Compression works",{
  mat = array(1:(4^4), dim=rep(4,4))
  
  cmat = compress.data(mat)

  cmat.brute = array(NA, dim=c(rep(2,3),4))
  for(t in 1:4){
    for(i in 1:2){
      for(j in 1:2){
        for(k in 1:2){
          cmat.brute[i,j,k,t] = mean(mat[(2*(i-1)+1):(2*i), (2*(j-1)+1):(2*j),
           (2*(k-1)+1):(2*k), t])
        }
      }
    }
  }

  expect_true(all(cmat.brute == cmat))
})

#test if I compress a matrix with even dimensions, the flip of the matrix
#  should give me the flip of the output
revengen <- function(){
  dim = rinteger(elements = {size = c(min = 1, max = 24)}, size = c(min = 3, max = 4))
  dim[1:3] = 2*dim[1:3]
  x = rdouble(size = c(min = prod(dim), max = prod(dim)))

  mat = array(x, dim=dim)
  mat.rev = array(rev(x), dim=dim)
}
test.out = ##INCOMPLETE
