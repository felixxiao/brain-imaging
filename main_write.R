files = grep('-criteria.RData', dir('ABIDE'), value = T)
files = paste0('ABIDE/', files)

criteria = c('Adjacent', 'Boundary',  'RatioCut',
             'CompParc',  'Balance', 'Jaggedness')
methods  = c('GenEC (3,1)', 'GenEC (6,1)', 'GenEC (6,4)',
             'Spectral', 'Spectral GenEC(6,4)', 'SymBMF', 'SymBMF GenEC(6,4)')
brains   = sapply(strsplit(files, '[[:punct:]]'), function(x) x[2])
crit = array(NA, dim = c(length(criteria), length(brains), length(methods)),
             dimnames = list(criteria, brains, methods))
for (i in 1:length(files))
{
  load(files[i])
  crit['Adjacent',  brains[i], ] = adjcent
  crit['Boundary',  brains[i], ] = boundry
  crit['RatioCut',  brains[i], ] = ratioct
  crit['CompParc',  brains[i], ] = connect
  crit[ 'Balance',  brains[i], ] = balance
  crit['Jaggedness',brains[i], ] = jaggedn
}

crit.mean = t(round(apply(crit[,brains[7:12],], c(1,3), mean), digits = 3))
write.table(crit.mean, file = 'writeup/figs/8_control.tex', quote = F,
            sep = ' & ', eol = ' \\\\\n')

methods = as.matrix(read.csv('writeup/figs/8_methods.csv'))
methods = plyr::mapvalues(methods, c('Y', 'N', 'Mostly'),
                          c('\\cellcolor{green!25}Yes',
                            '\\cellcolor{red!25}No',
                            '\\cellcolor{yellow!25}Mostly'))
write.table(methods, 'writeup/figs/8_methods.tex', quote = F, sep = ' & ',
            eol = ' \\\\\n', row.names = F)
