# Arguments:
#   edge_mat.csv  weights.csv  alpha  beta  500 400 300 ...

execfile('partition_edgecontract/source_contraction.py')

edge_mat = np.genfromtxt(sys.argv[1], dtype = int, delimiter = ',')
weights  = np.genfromtxt(sys.argv[2])

alpha = float(sys.argv[3])
beta  = float(sys.argv[4])

num_components = [int(k) for k in sys.argv[5:]]

cg = ContractibleGraph(edge_mat, weights)

# minimum priority is selected
def priority_func(cg, A):
  size = cg.get_size(A)
  weights = cg.get_links(A)
  bounds = cg.get_boundaries(A)
  priorities = {B: weights[B]**alpha * bounds[B] for B in weights}
  B = max(priorities, key = priorities.get)
  return (size ** beta / priorities[B], B)

write_partition_csv(cg, num_components, priority_func)
