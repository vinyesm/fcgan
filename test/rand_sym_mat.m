function [ M ] = rand_sym_mat( sz,rank )

spec=[abs(rand(1,rank)) zeros(1,sz-rank)];
[Q,R]=qr(randn(sz));
M=Q'*diag(spec)*Q;

end

