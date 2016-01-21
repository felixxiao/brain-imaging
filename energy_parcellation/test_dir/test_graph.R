#make a dummy dataset that we'll consider
#construction: it's a 5x5x6x10 dataset. The last
# 10 is for time.
#Take a look at the 5x5x6 box: (called X)
# the elements X[,,1], X[,,5], X[,,6],
#              X[1,,], X[5,,],
#              X[,1,], X[,5,] are all 0 (empty)
# the elements X[2:3, 2:3, 1:4] are drawn from one distribution
# the elements X[4, 4, 1:4] are drawn from another
# the elements X[2:3, 4, 1:4] another, and X[4, 2:3, 1:4]


