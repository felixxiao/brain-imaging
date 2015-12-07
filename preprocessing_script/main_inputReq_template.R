#!/usr/bin/Rscript

source("source_header.R")

parser = ArgumentParser()
parser$add_argument("-t", "--template", type="character", 
 help="File location of template .nii.gz file")
parser$add_argument("-o", "--output", type="character",
 help="File prefix on where to save the output")

args = parser$parse_args()

template = readNIfTI(args$template)
template = template@.Data

fileoutput = paste0(args$output, "_", Sys.Date(), ".RData")

res = extract.neighbors(template, tmpsave = fileoutput)
save(res, file = fileoutput)
