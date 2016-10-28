function [v,i]=max_k_chain(x,k)
% For a vector x in R^p computes the index of the contiguous block of 
% length k that has maximal norm.
% The block with maximal norm is [i,...,i+k-1] where i is returned
m=size(x,2);
if (m~=1),
    error('x is not a column vector');
end
z=conv(x.^2, ones(k,1),'valid');
[v,i]=max(z);
v=sqrt(v);

