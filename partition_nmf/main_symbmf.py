import numpy as np
import sys
from scipy import sparse

# A is either numpy.array or scipy.sparse.csc_matrix
def SymBMF_local(A, k, incomplete = False):
	n = A.shape[1]

	# initialize X : n by k
	X = map(lambda j: [0]*j + [1] + [0]*(k-j-1), np.random.randint(0, k, n))
	X = np.array(X)

	def cost(A, i, X):
		if type(A) is sparse.csc_matrix:
			a = A.getcol(i)
			idx = a.nonzero()
			return np.sum(np.abs(a[idx].T - X[idx[0],:]), 0)
		else:
			return np.sum(np.abs(A[:,i] - X.T), 1)

	i = 0
	j = 0
	while True:
		k = np.argmax(X[i,:])
		X[i,k] = 0

		k_best = np.argmin(cost(A, i, X))

		X[i,k_best] = 1

		if k != k_best:
			j = i
		i = (i + 1) % n
		if i == j:
			break

	return X


if len(sys.argv) != 1 + 3:
	print 'Arguments: A.csv  k  X.csv'
	sys.exit()

A = np.genfromtxt(sys.argv[1], delimiter = ',')
A = sparse.csc_matrix(A)

k = int(sys.argv[2])

X, obj = SymBMF_local(A, k)

np.savetxt(sys.argv[3], X, fmt = '%d', delimiter = ',')
