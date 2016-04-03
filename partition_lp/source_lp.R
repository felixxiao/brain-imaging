write.ampl_data.fractional = function(file, edges, k, n = NULL,
                                      assign.mat = T)
{
  if (is.null(n)) n = length(unique(as.vector(edges$edge.mat)))
  m = nrow(edges$edge.mat)

  sink(file)

  cat('param m := ', m, ';\n\n', sep = '')
  cat('param K := ', k, ';\n\n', sep = '')

  if (assign.mat)
  { 
    cat('param n := ', n, ';\n\n', sep = '')

    cat('param h :=\n')
    for (j in 1:m)
      cat(sprintf('%8d %8d\n', j, edges$edge.mat[j, 1]))
    cat(';\n\n')

    cat('param i :=\n')
    for (j in 1:m)
      cat(sprintf('%8d %8d\n', j, edges$edge.mat[j, 2]))
    cat(';\n\n')
  }

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
  out = sapply(rows, function(r) as.numeric(r[2:length(r)]))
  if (is.matrix(out)) return(t(out))
  else                return(out)
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

partition.fractional.approx = function(edges, k)
{
  write.ampl_data.fractional('data_fraction_approx.ampl', edges, k,
                             assign.mat = F)
  #output = system('ampl fractional_approx.ampl', intern = T)


}