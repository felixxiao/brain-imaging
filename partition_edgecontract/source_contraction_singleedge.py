import numpy as np
import heapq
import sys
from datetime import datetime
import csv
import json

class ContractibleGraph:
  # assumption : every pair of vertices have at most one edge
  # assumption : no vertex has an edge to itself
  # edge_mat   : numpy.array of shape (|E|, 2), no redundancies
  def __init__(self, edge_mat = None, weights = None):
    if edge_mat is None:
      self.vertices = {}
      self.N = 0
      self.edges = {}
      return
    self.vertices = {v: [v] for v in np.unique(edge_mat)}
    self.N = len(self.vertices)
    self.edges = {v: {} for v in self.vertices}

    for i in xrange(edge_mat.shape[0]):
      v, w = edge_mat[i,:]
      # assert v != w
      self.edges[v][w] = self.edges[w][v] = [weights[i], 1]

  @classmethod
  def read_file(cls, filename):
    f = open(filename, 'r')
    objects = json.load(f)
    f.close()
    cg = cls()
    cg.vertices = {int(C): V for C, V in objects['vertices'].items()}
    cg.N        = objects['N']
    cg.edges    = {int(A): {int(B): link for B, link in links.items()}
                   for A, links in objects['edges'].items()}
    return cg

  def contract_components(self, A, B):
    # assert A and B in self.vertices.keys()
    # assert A in self.edges[B].keys()
    # assert B in self.edges[A].keys()

    self.vertices[A].extend(self.vertices[B])
    self.vertices.pop(B)

    for B_adj in self.edges[B].keys():
      if B_adj in self.edges[A].keys():
        # assert A in self.edges[B_adj].keys()
        # assert self.edges[B_adj][A] is self.edges[A][B_adj]
        self.edges[A][B_adj][0] += self.edges[B][B_adj][0]
        self.edges[A][B_adj][1] += self.edges[B][B_adj][1]
      elif B_adj != A:
        # assert A not in self.edges[B_adj].keys()
        self.edges[A][B_adj] = self.edges[B_adj][A] = list(self.edges[B][B_adj])
      self.edges[B_adj].pop(B)

    self.edges.pop(B)

  def get_links(self, A):
    return {B: link[0] / link[1] for B, link in self.edges[A].items()}

#  def get_link_weight(self, A, B):
#    return self.edges[A][B][0] / self.edges[A][B][1]

  def get_size(self, A):
    return len(self.vertices[A])

  def get_boundaries(self, A):
    return {B: cg.edges[A][B][1] / min(cg.get_size(A), cg.get_size(B))
            for B in self.edges[A]}

  def get_vertex_components(self):
    return {v: k for k, V in self.vertices.items() for v in V}

  def get_vertices_sorted(self):
    V = [v for C in self.vertices.values() for v in C]
    return sorted(V)

  def save_file(self, filename = None):
    if filename is None:
      filename = 'cg' + str(len(self.vertices)) + '.json'
    f = open(filename, 'w')
    objects = {'vertices' : cg.vertices,
               'N'        : cg.N,
               'edges'    : cg.edges}
    json.dump(objects, f, separators = (',',':'))
    f.close()

# priority_func(cg, A) returns a tuple (priority, endpoint)
# minimum priority is selected
def partition_contractedge(num_components, cg, priority_func,
                           disp_all = False):
  k = len(cg.vertices)
  assert num_components > 0
  if k <= num_components:
    print 'CG already has target number of components'
    return
  print 'Computing priorities'

  priority = {A: priority_func(cg, A) for A in cg.vertices}
  pq = [(x[0], A) for A, x in priority.items()]
  heapq.heapify(pq)

  print 'Contracting graph'
  start_time = datetime.now()
  while k > num_components:
    val, A = heapq.heappop(pq)
    if A in cg.vertices:
      priority = priority_func(cg, A)
      if val != priority[0] and val > pq[0][0]: # PQ entry out-of-date
        heapq.heappush(pq, (priority[0], A))
      else:
        B = priority[1]
        cg.contract_components(A, B)
        k -= 1
        if k == 1: break

        priority = priority_func(cg, A)
        heapq.heappush(pq, (priority[0], A))
        
        if k % 10000 == 0:
          print 'Number of components:', k
          print '  len(pq)', len(pq)
          print '  Time:', datetime.now() - start_time
          cg.save_file()
          start_time = datetime.now()
  cg.save_file()

# num_components must be strictly decreasing
def write_partition_csv(cg, num_components, priority,
                        filename = 'partition.csv'):
  f = open(filename, 'w')
  writer = csv.writer(f)
  writer.writerow(num_components)
  V = cg.get_vertices_sorted()
  writer.writerow(V)
  for k in num_components:
    partition_contractedge(k, cg, priority)
    partition = cg.get_vertex_components()
    partition = [partition[v] for v in V]
    writer.writerow(partition)
  f.close()

#"""
if __name__ == '__main__':
  edge_mat = np.genfromtxt('edge_mat.csv', dtype = int, delimiter = ',')
  weights  = np.genfromtxt('weights.csv')
  
  num_components = [int(k) for k in sys.argv[1:]]

  cg = ContractibleGraph(edge_mat, weights)

  f = lambda cg, A: {B: cg.get_size(B) * (1 \
                      - cg.get_link_weight(A, B) \
                      - cg.get_boundary(A, B)) for B in cg.edges[A]}

  # minimum priority is selected
  def priority_func(cg, A):
    size = cg.get_size(A)
    weights = cg.get_links(A)
    bounds = cg.get_boundaries(A)
    priorities = {B: weights[B]**6 * bounds[B] for B in weights}
    B = max(priorities, key = priorities.get)
    return ((size + 10) / priorities[B], B)

  write_partition_csv(cg, num_components, priority_func)

"""
edge_mat = np.array([[1, 2],
                     [1, 3],
                     [2, 4],
                     [3, 4]])
weights = [0.7, 0.8, 0.5, 0.9]

edge_mat = np.array([[1, 2],
                     [1, 3],
                     [1, 4],
                     [2, 3],
                     [3, 5],
                     [3, 6],
                     [4, 5],
                     [5, 6]])
weights = [0.9, 0.8, 0.69, 0.5, 0.6, 0.4, 0.7, 0.6]

cg = ContractibleGraph(edge_mat, weights)

"""