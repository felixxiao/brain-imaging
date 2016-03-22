import numpy as np
import heapq
import sys

class ContractibleGraph:
  # assumption : every pair of vertices have at most one edge
  # assumption : no vertex has an edge to itself
  # edge_mat   : numpy.array of shape (|E|, 2), no redundancies
  def __init__(self, edge_mat, weights):
    self.vertices = {v: [v] for v in np.unique(edge_mat)}
    self.N = len(self.vertices)
    self.edges = {v: {} for v in self.vertices}
    
    for i in xrange(edge_mat.shape[0]):
      v, w = edge_mat[i,:]
      assert v != w
      self.edges[v][w] = self.edges[w][v] = [weights[i]]

  def contract_components(self, A, B):
    assert A and B in self.vertices.keys()
    assert A in self.edges[B].keys()
    assert B in self.edges[A].keys()

    self.vertices[A].extend(self.vertices[B])
    self.vertices.pop(B)

    for B_adj in self.edges[B].keys():
      if B_adj in self.edges[A].keys():
        assert A in self.edges[B_adj].keys()
        assert self.edges[B_adj][A] is self.edges[A][B_adj]
        self.edges[A][B_adj].extend(self.edges[B][B_adj])
      elif B_adj != A:
        assert A not in self.edges[B_adj].keys()
        self.edges[A][B_adj] = self.edges[B_adj][A] = list(self.edges[B][B_adj])
      self.edges[B_adj].pop(B)

    self.edges.pop(B)

  def get_links(self, A):
    return {B: np.mean(weights) for B, weights in self.edges[A].items()}

  def get_vertex_components(self):
    return {v: k for k, V in self.vertices.items() for v in V}

def partition_contractedge(num_components, cg):
  print 'Computing priorities'
  priority = {A: (len(cg.vertices[A]) - max(links),
                  max(links, key = links.get))
              for A, links in {A: cg.get_links(A) for A in cg.vertices}.items()}
  pq = [(x[0], A) for A, x in priority.items()]
  heapq.heapify(pq)

  print 'Contracting graph'
  k = cg.N
  while k > num_components:
    val, A = heapq.heappop(pq)
    if A in priority:
      if priority[A][0] != val:
        heapq.heappush(pq, (priority[A][0], A))
      else:
        B = priority[A][1]
        assert B in priority
        cg.contract_components(A, B)
        k -= 1

        for C in cg.get_links(A):
          assert C in priority
          links = cg.get_links(C)
          priority[C] = (len(cg.vertices[C]) - max(links),
                         max(links, key = links.get))
        links = cg.get_links(A)
        priority[A] = (len(cg.vertices[A]) - max(links),
                       max(links, key = links.get))
        heapq.heappush(pq, (priority[A][0], A))
        priority.pop(B)
        if k % 10000 == 0: print k,
  return cg

if __name__ == '__main__':
  # Arguments: edge_mat_file.csv  weights_file.csv  num_components
  edge_mat = np.genfromtxt(sys.argv[1], dtype = int, delimiter = ',')
  weights  = np.genfromtxt(sys.argv[2])
  print 'Initializing CG'
  cg = ContractibleGraph(edge_mat, weights)

  cg_new = partition_contractedge(int(sys.argv[3]), cg)

  assert cg_new is cg

"""
edge_mat = np.array([[1, 2],
                     [1, 3],
                     [2, 4],
                     [3, 4]])
weights = [0.7, 0.8, 0.5, 0.9]
cg = ContractibleGraph(edge_mat, weights)

cg.get_links(1)
"""