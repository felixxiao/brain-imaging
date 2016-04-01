set V;

param A {V,V};         # adjacency matrix
param n;
param k;
param MIN_SIZE := 1;
param MAX_SIZE := n;

var Z {V,V} >= 0;      # partition matrix
var X {V,1..k} binary; # assignment matrix

maximize within:
	sum {i in V, j in V} Z[i,j] * A[i,j];

subject to exclusive {i in V}:
  sum {h in 1..k} X[i,h] = 1;

subject to min_size {h in 1..k}:
  sum {i in V} X[i,h] >= MIN_SIZE;

subject to max_size {h in 1..k}:
  sum {i in V} X[i,h] <= MAX_SIZE;

subject to partit_assign_1 {i in V, j in V, h in 1..k}:
  Z[i,j] <= 1 + X[i,k] - X[j,k];

subject to partit_assign_2 {i in V, j in V, h in 1..k}:
  Z[i,j] <= 1 - X[i,k] + X[j,k];

data;

data data.ampl;

option cplex_options 'outlev=0';
option solver_msg 0;

solve;

#for {i in V}
#{
#	printf {j in V} "%f,", Z[i,j];
#	printf "\n";
#};

for {i in V}
{
  printf {h in 1..k} "%f,", X[i,h];
  printf "\n";
}
