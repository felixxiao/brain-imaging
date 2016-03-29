# Arguments
#   cg.json  alpha  beta  500 400 300 ...

execfile('source_contraction.py')

cg = ContractibleGraph.read_file(sys.argv[1])

alpha = float(sys.argv[2])
beta  = float(sys.argv[3])

num_components = [int(k) for k in sys.argv[4:]]

# minimum priority is selected
def priority_func(cg, A):
  size = cg.get_size(A)
  weights = cg.get_links(A)
  bounds = cg.get_boundaries(A)
  priorities = {B: weights[B]**alpha * bounds[B] for B in weights}
  B = max(priorities, key = priorities.get)
  return (size ** beta / priorities[B], B)

write_partition_csv(cg, num_components, priority_func)
