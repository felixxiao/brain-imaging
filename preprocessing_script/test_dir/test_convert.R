library(testthat)
setwd("..")
source("header_preprocess.R")

test_that("Conversion works both ways", {
  set.seed(10)
  box = array(0, dim=c(5,5,5))
  idx = sample(1:(5^3),20)
  box[idx] = 1
  dimen = dim(box)

  mask = which(box!=0)
  expect_true(all(sort(idx)==sort(mask)))

  for(i in 1:length(mask)){
    loc = convert.2Dto3Dloc(mask[i], dimen)
    expect_true(box[loc[1], loc[2], loc[3]]==1)
  }

  for(i in 1:dimen[1]){
    for(j in 1:dimen[2]){
      for(k in 1:dimen[3]){
        if(box[i,j,k]==1){
          idx = convert.3Dto2Dloc(c(i,j,k), dimen)
          expect_true(idx %in% mask)
        } 
      }
    }
  }
})

test_that("Extract neighbors works", {
  box = array(0, dim=c(5,5,5))
  box[c(2:3), c(1:3), c(4:5)] = 1
  dimen = dim(box)

  res = extract.neighbors(box)
  expect_true(length(res$mask) == (2*3*2))
  expect_true(length(res$mask) == length(res$neighbor.list))

  for(i in 1:length(res$mask)){
    idx = res$mask[i]
    loc = convert.2Dto3Dloc(idx, dimen)
    neigh = res$neighbor.list[[i]]

    for(j in 1:length(neigh)){
      loc2 = convert.2Dto3Dloc(res$mask[neigh[j]], dimen)
      expect_true(box[loc2[1], loc2[2], loc2[3]] == 1)
      expect_true(all(abs(loc2-loc)<=1))
    }
  }
})
