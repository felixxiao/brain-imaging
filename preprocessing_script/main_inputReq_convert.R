#!/usr/bin/Rscript

rm(list=ls())
source("header_preprocess.R")

parser = ArgumentParser()
parser$add_argument("-i", "--input", type="character",
 help="File location of input fMRI data in .nii.gz format")
parser$add_argument("-o", "--output", type="character",
 help="File prefix on where to save the output")
parser$add_argument("-t", "--template", type="character",
 help="File output from main_inputReq_template.R")
parser$add_argument("-d", "--details", type="character",
 help="Addition details to append to output list as a string")

args = parser$parse_args()
save(args,file="tmp.RData")


dat = readNIfTI(args$input)
dat = dat@.Data

load(args$template) #loads a file called "template"
mask = template$mask

fileoutput = paste0(args$output, "_", Sys.Date(), ".RData")

#save output as a sparse matrix
res = extract.data(dat, mask)
#res = Matrix(res, sparse=TRUE)
res = list(dat = res, details = args$details)

dat = res
save(dat, file = fileoutput)
