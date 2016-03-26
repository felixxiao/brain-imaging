param m;
param K;
param n;

param h {1..m};
param i {1..m};

param a {1..m};

var   w {1..m, 1..K} >= 0;
var   z {1..m, 1..K};
var   x {1..n, 1..K} binary;
var   y {1..m};

maximize mean_cut:
  sum {j in 1..m, k in 1..K} a[j] * w[j,k];


subject to quad_line_1 {j in 1..m, k in 1..K}:
  1 + z[j,k] >= x[h[j],k] + x[i[j],k];

subject to quad_line_2 {j in 1..m, k in 1..K}:
  z[j,k] <= x[h[j],k];

subject to quad_line_3 {j in 1..m, k in 1..K}:
  z[j,k] <= x[i[j],k];


subject to min_edges {k in 1..K}:
  sum {j in 1..m} z[j,k] >= 1;

subject to assign {l in 1..n}:
  sum {k in 1..K} x[l,k] == 1;

subject to min_size {k in 1..K}:
  sum {l in 1..n} x[l,k] >= 1;


subject to frac {k in 1..K}:
  sum {j in 1..m} w[j,k] == 1;


subject to bin_line_1 {j in 1..m, k in 1..K}:
  y[j] - w[j,k] <= 1 - z[j,k];

subject to bin_line_2 {j in 1..m, k in 1..K}:
  w[j,k] <= y[j];

subject to bin_line_3 {j in 1..m, k in 1..K}:
  w[j,k] <= z[j,k];

subject to bin_line_4 {j in 1..m, k in 1..K}:
  w[j,k] >= 0;


data;

solve;
display x;
