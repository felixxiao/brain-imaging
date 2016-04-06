import numpy as np
import heapq
import sys
from datetime import datetime
import csv
import json

class ContractibleGraph:
  # assumption : every pair of vertices have at most one edge
  # assumption : no vertex has an edge to itself
  # edge_mat   : numpy.array of shape (|E|, 2), redundancies added
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
      if v == w: continue
      if w not in self.edges[v].keys():
        assert v not in self.edges[w]
        self.edges[v][w] = self.edges[w][v] = [weights[i], 1]
      else:
        assert v in self.edges[w]
        self.edges[w][v][0] += weights[i]
        self.edges[w][v][1] += 1
      assert self.edges[w][v] is self.edges[v][w]

  @classmethod
  def read_file_json(cls, filename):
    f = open(filename, 'r')
    objects = json.load(f)
    f.close()
    cg = cls()
    cg.vertices = {int(C): V for C, V in objects['vertices'].items()}
    cg.N        = objects['N']
    cg.edges    = {int(A): {int(B): link for B, link in links.items()}
                   for A, links in objects['edges'].items()}
    for A in cg.edges:
#      print A
      assert A in cg.vertices
      assert A not in cg.edges[A]
      for B in cg.edges[A]:
        assert B in cg.vertices
        assert cg.edges[B][A] == cg.edges[A][B]
        if cg.edges[A][B] is not cg.edges[B][A]:
          cg.edges[A][B] = cg.edges[B][A]
    return cg

  def assert_valid(self):
    for A in self.edges:
      if A not in self.vertices:
        raise Exception(str(A) + ' not in vertices')
      for B in self.edges[A]:
        if B not in self.vertices:
          raise Exception(str(B) + ' not in vertices')
        if self.edges[B][A] is not self.edges[A][B]:
          raise Exception('edges between ' + str(A) + ' and ' + str(B) + ' unlinked')
        if A not in self.edges[B]:
          raise Exception(str(A) + ' not in edges[' + str(B) + ']')

  @classmethod
  def read_files(cls, edges_csv, vertices_json):
    edges = np.genfromtxt(edges_csv, delimiter = ',')
    edge_mat = edges[:,0:2].astype(int)
    weights  = edges[:,2]
    cg = cls(edge_mat, weights)
    f = open(vertices_json)
    cg.vertices = json.load(f)
    f.close()
    return cg

  def contract_components(self, A, B):
#    assert A and B in self.vertices
#    assert A in self.edges[B]
#    assert B in self.edges[A]
#    assert A != B
#    self.assert_valid()

    self.vertices[A].extend(self.vertices[B])

    for B_adj in self.edges[B].keys():
#      assert self.edges[B_adj][B] is self.edges[B][B_adj]
#      assert B_adj != B
      if B_adj in self.edges[A]:
#        assert A in self.edges[B_adj]
#        assert self.edges[B_adj][A] is self.edges[A][B_adj]
        self.edges[A][B_adj][0] += self.edges[B][B_adj][0]
        self.edges[A][B_adj][1] += self.edges[B][B_adj][1]
      elif B_adj != A: # B_adj not in self.edges[A]
        self.edges[A][B_adj] = self.edges[B_adj][A] = list(self.edges[B][B_adj])
      self.edges[B_adj].pop(B)
      self.edges[B].pop(B_adj)
#      self.assert_valid()

    self.vertices.pop(B)
    self.edges.pop(B)

  def get_links(self, A):
    return {B: link[0] / link[1] for B, link in self.edges[A].items()}

#  def get_link_weight(self, A, B):
#    return self.edges[A][B][0] / self.edges[A][B][1]

  def get_size(self, A):
    return len(self.vertices[A])

  def get_boundaries(self, A):
    return {B: float(cg.edges[A][B][1]) / min(cg.get_size(A), cg.get_size(B))
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
          start_time = datetime.now()

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
    cg.save_file()
    partition = cg.get_vertex_components()
    partition = [partition[v] for v in V]
    writer.writerow(partition)
  f.close()

