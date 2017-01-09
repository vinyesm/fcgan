function [uBest,kBest,allVal] = lmo_spsd(A,param)

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
        param.k=k;
        param.q=k;
        [u,lambda] = TPI(B,param);
        lambda=(u'*A*u)/param.cardfun(param.k);
        allVal(k)=lambda;
        if lambdaBest < lambda
            uBest = u;
            lambdaBest = lambda;
            kBest=k;
        end
    end
end

end