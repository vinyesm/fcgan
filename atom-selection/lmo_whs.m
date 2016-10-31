function [maxval,new_atom]=lmo_whs(s,param)
if size(s,2)>1,
    error('s is not a column vector');
end
p=size(s,1);
groups=param.group_mat;
sgroups=bsxfun(@times,s,groups);
sgroups=sqrt(sum(sgroups.^2,1));

[maxval,i]=max(sgroups);
J=find(groups(:,i));
new_i=(J)';
new_j=ones(length(J),1);
new_vals=s(new_i);
new_vals=new_vals./norm(new_vals);

new_atom=sparse(new_i,new_j,new_vals,p,1);
if full(new_atom)'*s<=0,
    error('atom wrong');
end
