function [uBest vBest lambdaBest] = MultiTPI(A,param)

lambdaBest = -inf;
uBest = randn(size(A,1),1); uBest = projectL0(uBest,param.k); uBest = uBest/norm(uBest);
vBest = randn(size(A,2),1); vBest = projectL0(vBest,param.q); vBest = vBest/norm(vBest);
%         fprintf('length %d\n',length(find(vBest)))
  
for i=1:param.stPtPowerIter
    [u v] = truncatedPowerSeq(A,param);
%    fprintf('length %d\n',length(find(v)))
    lambda = (u'*A*v);
    if lambdaBest < lambda
        uBest = u;
        vBest = v;
        lambdaBest = lambda;
        
    end    
end

end