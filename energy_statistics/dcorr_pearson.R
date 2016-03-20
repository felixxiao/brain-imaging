library(MASS)
library(energy)

N = 2000
dist.cor = c()

for (rho in seq(-1, 1, by = 0.1))
{
  cat(rho, ' ')
  Sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  X = mvrnorm(N, c(0, 0), Sigma)
  
  #plot(X)
  dist.cor = c(dist.cor, dcor(X[,1], X[,2]))
}
plot(dist.cor)
