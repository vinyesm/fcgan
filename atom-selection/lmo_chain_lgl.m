function [maxval,new_atom]=lmo_chain_lgl(s,param)
if size(s,2)>1,
    error('s is not a column vector');
end
p=size(s,1);
[maxval,i]=max_k_chain(s,param.K);
new_i=(i:i+param.K-1)';
new_j=ones(param.K,1);
new_vals=s(new_i);
new_vals=new_vals./norm(new_vals);

new_atom=sparse(new_i,new_j,new_vals,p,1);
if full(new_atom)'*s<=0,
    error('atom wrong');
end
