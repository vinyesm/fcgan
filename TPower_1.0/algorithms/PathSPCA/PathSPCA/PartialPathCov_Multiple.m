function [ V ]=PartialPathCov_Multiple(S, sparsity)
% Given covariance matrix, compute full sparse PCA path
% Input: 
%   S: covariance matrix with decreasing diagonal
% Output: 
%   vars: vector of variances for each target cardinality
%   rhobreaks: the coresponding rho penalties
%   res: each column contains the corresponding subset of variables

[d,ix]=sort(diag(S),'descend');
S_sort = S(ix,ix);

ds=diag(S_sort);
% if any(ds(1:end-1)-ds(2:end)<0)
%     disp('Error in FullPathCov input: diagonal of input matrix should be decreasing');
%     isopt=0;rho=NaN;return;
% end



% Loop through variables
V = [];
level = length(sparsity);

for (l = 1:level)
    sparsity_cur = sparsity(l);
    
    n=size(S_sort,1);A=chol(S_sort);
    subset=[1];subres=[subset';zeros(n-length(subset),1)];
    res=[];rhobreaks=[sum(A(:,1).^2)];sol=[];vars=[];

    for i=1:sparsity_cur
        % Compute solution at current subset
        [v,mv]=maxeig(S_sort(subset,subset));
        vsol=zeros(n,1);vsol(subset)=v;
        sol=[sol,vsol];vars=[vars,mv];
        % Compute x at current subset
        x=A(:,subset)*v;x=x/norm(x);
        res=[res,[subset';zeros(n-length(subset),1)]];
        % Compute next rho breakpoint
        set=1:n;set(subset)=[];
        vals=(x'*A(:,set)).^2;
        [rhomax,vpos]=max(vals);
        rhobreaks=[rhobreaks;rhomax];
        subset=[subset,set(vpos)];
    end
    vsol_ori = vsol;
    vsol_ori(ix) = vsol;
    V = [V , vsol_ori];
    %deflation
    if(level>1)
        P = eye(n)- vsol_ori*vsol_ori';
        S = P*S*P;
        S = S+eye(size(S))*1e-8;
        
        [d,ix]=sort(diag(S),'descend');
        S_sort=S(ix,ix);

    end
end