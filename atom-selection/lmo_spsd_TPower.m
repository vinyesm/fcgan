function [uBest,kBest,allVal] = lmo_spsd_TPower(A,param)

B=0.5*(A+A');
% emin=eigs(A,1,'sa');
% if emin>0
%     emin=0;
% end
% B=A-1.1*emin*eye(size(A,1));
lambdaBest= -inf;
kBest=0;
allVal=zeros(size(A,1),1);

for k=1:size(A,1)
    if (param.cardfun(k) ~= inf)
        options.verbose=0;
        options.optTol=1e-8;
        options.maxIter=1000;
        options.cardinality_vec=k;
        [u,lambda] = TPower_SPCA(A, options);
        lambda=lambda/param.cardfun(k);
        allVal(k)=lambda;
        if lambdaBest < lambda
            uBest = u;
            lambdaBest = lambda;
            kBest=k;
        end
    end
end

end