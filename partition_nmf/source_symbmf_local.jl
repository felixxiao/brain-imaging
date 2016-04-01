import StatsBase.sample
using DataFrames
using Convex

function partition_opt_convex(L, k)
	n = size(L, 1)
	Z = Variable(n, n, Positive())
	problem = minimize(trace(Z * L), [
	  diag(Z) == 1,
	  sum(Z, 1) == n / k,
	  Z <= 1,
	  Z == Z'
	])

	solve!(problem)
	return Z.value
end

function factor_01(Z, k; initial = 0, ITER = 4)
	n = size(Z, 1)
	if initial == 0
		X = falses(n,k)
		s = StatsBase.sample(1:k, n, replace = true)
		for i = 1:n
	  	X[i, s[i]] = true
		end
	else
		@assert typeof(initial) <: Array{Bool}
		X = copy(initial)
	end

	for i = 1:n
	  Z[i,i] = 0
	end

	r = ones(k)
	for l = 1:ITER
	  println("Iter ", l, " cost ", vecnorm(Z - X * X'))
	  for i = 1:n
	    X[i,:] = false
	    for j = 1:k
	      r[j] = sum(abs(Z[:,i] - X[:,j]))
	    end
	    j = findmin(r)[2]
	    X[i,j] = true
	  end
	end
	println()

	for i = 1:n
	  Z[i,i] = 1
	end

	return (X, vecnorm(Z - X * X'))
end

if length(ARGS) < 3
	println("Arguments: Z.csv  k  iter")
end

Z = readtable(ARGS[1], header = false)
Z = Array(Z[:, 1:nrow(Z)])

# number of components
k = parse(ARGS[2])

X, c = factor_01(Z, k, ITER = parse(ARGS[3]))

writetable("X.csv", DataFrame(X*1), header = false)
