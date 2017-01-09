function [vars_adj_rate] = vars_adj(u, A)

B = u'*A*u;
[U D V] = svd(B);
Z = U*D.^0.5;

[Q R] = qr(Z);
r = diag(R);
var_adj = norm(r)^2;
vars_adj_rate = var_adj / trace(A);