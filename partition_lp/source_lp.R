write.ampl_data.fractional = function(file, edges, k, n = NULL)
{
  if (is.null(n)) n = max(edges$edge.mat)
  m = nrow(edges$edge.mat)

  sink(file)

  cat('param m := ', m, ';\n\n', sep = '')
  cat('param K := ', k, ';\n\n', sep = '')
  cat('param n := ', n, ';\n\n', sep = '')

  cat('param h :=\n')
  for (j in 1:m)
    cat(sprintf('%8d %8d\n', j, edges$edge.mat[j, 1]))
  cat(';\n\n')

  cat('param i :=\n')
  for (j in 1:m)
    cat(sprintf('%8d %8d\n', j, edges$edge.mat[j, 2]))
  cat(';\n\n')

  cat('param a :=\n')
  for (j in 1:m)
    cat(sprintf('%8d %f\n', j, edges$energy.vec[j]))
  cat(';\n')

  sink()
}

read.ampl_output = function(file, start_char = ':', end_char = ';')
{
  lines = readLines(file)
  start_line = which(sapply(lines, substr, start = 1, stop = 1) ==
                     start_char) + 1
  end_line   = which(sapply(lines, substr, start = 1, stop = 1) ==
                     end_char) - 1
  rows = strsplit(lines[start_line:end_line], split = '\\s+')
  t(sapply(rows, function(r) as.numeric(r[2:length(r)])))
}

write.ampl_data = function(file, mat, k, laplacian = FALSE,
                           min_size = NULL, max_size = NULL)
{
  n = nrow(mat)

  sink(file)

  cat('param n := ', n, ';\n\n', sep = '')
  cat('param k := ', k, ';\n\n', sep = '')

  if (! is.null(min_size))
    cat('param MIN_SIZE := ', min_size, ';\n\n', sep = '')
  if (! is.null(max_size))
    cat('param MAX_SIZE := ', max_size, ';\n\n', sep = '')

  cat('set V := ')
  cat(1:n)
  cat(';\n\n')

  if (laplacian) cat('param L : ')
  else           cat('param A : ')
  cat(sprintf('%8d', 1:n))
  cat(' :=\n')
  for (i in 1:n)
  {
    cat(sprintf('%9d ', i))
    cat(sprintf('%.6f', mat[i,]))
    cat('\n')
  }
  cat(';\n\n')

  sink()
}

# needs work
compute.factor_01.julia = function(csv.file, num.components, iter = 5)
{
  julia.file = 'partition_opt/source_factor.jl'
  shell(paste('julia', julia.file, csv.file, num.components, iter))

  X = read.table('partition_opt/X.csv', sep = ',')
  parcel = rep(NA, times = nrow(X))
  for (j in 1:ncol(X)) parcel[which(X[,j])] = j

  as.factor(parcel)
}

compute.factor_qp = function(Z, X_init, lambda, ITER = 10)
{
#  Z = as.matrix(read.csv(csv.file, header = F))
  assert_that(isSymmetric(Z))
  assert_that(class(lambda) %in% c('function', 'numeric'))
  n = nrow(Z)
  X = X_init
  k = ncol(X)
  cost = rep(NA, times = ITER)

  if (class(lambda) == 'function')
    lambda = sapply(1:ITER, lambda)
  else if (length(lambda) == 1)
    lambda = rep(lambda, times = ITER)
  else
    assert_that(length(lambda) == ITER)

  for (t in 1:ITER)
  {
    cat(t, ' ')
    # \| Z - X X^T \|_F^2
    cost[t] = norm(Z - X %*% t(X), 'F')

    # min. 1/2 x^T Q x - p^T x
    p = as.vector(crossprod(X, Z + lambda[t] * diag(n)))
    Q = kronecker(diag(n), crossprod(X)) + lambda[t] * diag(n * k)

    # s.t. A^T x >= b
    A = - kronecker(diag(n), matrix(1, nrow = k))
    A = cbind(A, diag(n*k))
    b = c(rep(-1, times = n), rep(0, times = n*k))

    opt = solve.QP(Q, p, A, b, n)
    X = t(matrix(opt$solution, nrow = k))
  }
  X
}
