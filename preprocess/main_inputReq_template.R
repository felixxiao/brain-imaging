#!/usr/bin/Rscript

rm(list=ls())
source("header_preprocess.R")

parser = ArgumentParser()
parser$add_argument("-t", "--template", type="character", 
 help="File location of template .nii.gz file")
parser$add_argument("-o", "--output", type="character",
 help="File prefix on where to save the output")
parser$add_argument("-p", "--pattern", type="double", default=27,
 help="Number of voxels in a neighboring group. Either 7 or 27")
parser$add_argument("-d", "--details", type="character",
 help="Addition details to append to output list as a string")

args = parser$parse_args()

assert_that(args$pattern == 7 | args$pattern == 27)
if(args$pattern == 7){
  args$pattern = cross.enumerate()
} else {
  args$pattern = cube.enumerate()
}

template = readNIfTI(args$template)
template = template@.Data

fileoutput = paste0(args$output, "_", Sys.Date(), ".RData")

res= extract.neighbors(template, pattern = args$pattern,
 tmpsave = fileoutput)
res$details = args$details

template = res
save(template, file = fileoutput)
