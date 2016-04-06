execfile('source_symnmf.py')

if len(sys.argv) != 1 + 2:
	print 'Arguments: edges.txt  k'
	sys.exit()

X = partition_SymBMF_local(sys.argv[1], int(sys.argv[2]),
	                         incomplete = True)
np.savetxt('X_symbmf.csv', X, fmt = '%d', delimiter = ',')
