Felix's Senior Thesis (2015-16)
===============================
Code repository by Felix Xiao and Kevin Lin, under guidance of Han Liu.

#Folder Directory
##preprocessing\_script
(By Kevin) Code to convert a preprocessed FMRI data into a 2D matrix and an adjacency list

##common\_func
(By Felix and Kevin) A broad set of basic and generic functions (plotting, compressing data, etc.)
These functions are not specific for partitioning the brain, but are more for data processing and visualization.

##energy\_parcellation
(By Felix) For various functions to do parcellation using energy statistics

#Installation (Dependencies)
To run all the code in this repository, you will need the following R packages:
oro.nifti, assertthat, testtthat, argparse, igraph

In addition, you will need the following R packages via Git:
quickcheck

To install the latter set of packages, you can run the following lines
```
install.packages("devtools")
library("devtools")
install_github("RevolutionAnalytics/quickcheck@3.5.0", subdir = "pkg")
```

#How to Use
(To come)
