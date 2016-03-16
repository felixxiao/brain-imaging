set V;

param L {V,V};
param n;
param k;

var Z {V,V} >= 0;


minimize cost:
	sum {i in V, j in V} Z[i,j] * L[i,j];

subject to identity {i in V}:
	Z[i,i] = 1;

subject to size {i in V}:
	sum {j in V} Z[i,j] = n / k;

subject to upper {i in V, j in V}:
	Z[i,j] <= 1;

subject to symmetric {i in V, j in V}:
	Z[i,j] = Z[j,i];

data;

data data.ampl;

option cplex_options 'outlev=0';
option solver_msg 0;

solve;

for {i in V}
{
	printf {j in V} "%f,", Z[i,j];
	printf "\n";
};
