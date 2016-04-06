import numpy as np
from scipy import sparse
from numpy import matlib
from collections import defaultdict
import sys

def read_edges(filename):
  A = np.loadtxt(filename)
  w = A[:,2]
  i = A[:,0].astype(int) - 1
  j = A[:,1].astype(int) - 1
  w = np.concatenate((w, w))
  new_i = np.concatenate((i, j))
  j = np.concatenate((j, i))
  n = np.max(new_i)
  return sparse.csc_matrix((w, (new_i, j)))
  
def nonnegative_least_squares(CTC, CTB):
  assert type(CTC) and type(CTB) is np.matrix
  n = CTB.shape[1]
  k = CTC.shape[0]
  assert n > k
  assert k == CTC.shape[1] == CTB.shape[0]

  F = matlib.zeros((k, n), dtype = bool)
  G = matlib.ones((k, n), dtype = bool)
  I = range(n)
  X = matlib.zeros((k, n))
  Y = matlib.zeros((k, n))

  alpha = np.ones(n, dtype = int)
  beta  = np.ones(n, dtype = int) * (k + 1)

  def update(X, Y, F, G, I):
    group = defaultdict(list)
    for i in I:
      key = ' '.join(map(str, np.where(F[:,i])[0].tolist()[0]))
      group[key].append(i)
    for f, cols in group.items():
      if len(f) == 0:
        X[:,cols] = 0
        Y[:,cols] = - CTB[:,cols]
      else:
        f = [int(j) for j in f.split(' ')]
        X[np.ix_(f,cols)] = np.linalg.solve(CTC[np.ix_(f,f)],
                                            CTB[np.ix_(f,cols)])
        Y[np.ix_(f,cols)] = 0
        if len(f) < X.shape[0]:
          g = np.setdiff1d(np.arange(X.shape[0]), f)
          X[np.ix_(g,cols)] = 0
          Y[np.ix_(g,cols)] = CTC[np.ix_(g,f)] * X[np.ix_(f,cols)] \
                            - CTB[np.ix_(g,cols)]
    assert np.all(X[G] == 0) and np.all(Y[F] == 0)

  while True:
    # update feasible columns of X and Y
    update(X, Y, F, G, I)
    V = np.logical_or(X < 0, Y < 0)
    I = np.where(np.sum(V, axis = 0) > 0)[1].tolist()[0]
    if len(I) == 0:
      break
    for i in I:
      if np.sum(V[:,i]) < beta[i]:
        beta[i] = np.sum(V[:,i])
        alpha[i] = 3
      elif alpha[i] >= 1:
        alpha[i] -= 1
      else:
        j = np.max(np.where(V[:,i]))
        V[:,i] = np.zeros((k, 1), dtype = bool)
        V[j,i] = True
    F = np.logical_or(np.logical_xor(F, V), np.logical_and(V, G))
    G = np.logical_not(F)

  assert np.allclose(Y, CTC * X - CTB)
  assert np.all(Y >= 0)
  assert np.all(X >= 0)
  assert np.all(np.multiply(X, Y) == 0)

  return X

def partition_SymNMF(file_edges, k, tol = 1e-4,
                     alpha = 1, growth = .01):
  A = read_edges(file_edges)
  n = A.shape[0]

  HT = np.matrix(np.random.rand(k, n))

  t = 1
  while True:
    print t,
    WT = HT

    CTC = WT * WT.T + alpha * np.identity(k)
    CTB = WT * A    + alpha * WT

    HT = nonnegative_least_squares(CTC, CTB)
    alpha *= growth
    t += 1
    if np.max(np.abs(WT - HT)) <= tol:
      break

  HT.T

def partition_SymBMF_local(file_edges, k,
                           sparse = True, incomplete = False):
  if sparse:
    A = read_edges(file_edges)
  else:
    A = np.genfromtxt(file_edges, delimiter = ',')
  n = A.shape[1]

  # initialize X : n by k
  X = map(lambda j: [0]*j + [1] + [0]*(k-j-1), np.random.randint(0, k, n))
  X = np.array(X)

  if sparse:
    if incomplete:
      def cost(A, i, X):
        a = A.getcol(i)
        idx = a.nonzero()
        return np.sum(np.abs(a[idx].T - X[idx[0],:]), 0)        
    else:
      pass
  else:
    def cost(A, i, X):
      return np.sum(np.abs(A[:,i] - X.T), 1)

  i = 0
  j = 0
  while True:
    if i == 0: print '*',
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
