function [uBest,lambdaBest] = TPI(A,param)

lambdaBest = -inf;
uBest = randn(size(A,1),1); uBest = projectL0(uBest,param.k); uBest = uBest/norm(uBest);


for i=1:param.stPtPowerIter
    u = conGradU02(A,param.powerIter,param.k);
    lambda=(u'*A*u)/param.cardfun(param.k);
    if lambdaBest < lambda
        uBest = u;
        lambdaBest = lambda;
    end
end


end