execfile('source_symnmf.py')

if len(sys.argv) != 1 + 2:
  print 'Arguments: edges.txt  k'
  sys.exit()

X = partition_symnmf(sys.argv[1], int(sys.argv[2]))
np.savetxt('X_symnmf.csv', X, delimiter = ',')
