function [new_i,new_val,maxval]=get_new_atom_spca(Z,ActiveSet,param,inputData )

%
p=size(Z,1);
H = gradient(Z,inputData,param);
maxval=-inf;

for i=1:length(ActiveSet.I)
    %eigenvector associated to largest real eigeinvalue
    S=ActiveSet.I{i};
    Hs=H(S,S);
    Hs=0.5*(Hs+Hs');
    [v,d]=eigs(-Hs,1,'la');
    v=real(v);
    d=real(d);
    if d>maxval
        maxval=d;
        new_i=S;
        new_val=v;
%         maxatom=zeros(p,1);
%         maxatom(S)=v;
    end
end
% fprintf('\n%f',maxval);
if maxval==0
    display('Largest eigenvalue is negative or zero\n');
end


% anew=sparse(maxatom);


if ~isreal(v)
    error('selected new atom is not real\n')
end


end
