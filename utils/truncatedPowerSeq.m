function [u v] = truncatedPowerSeq(A,param)

debug=0;

v = randn(size(A,2),1);v = v/norm(v);
if param.PSD
    u = v;
else
    u = randn(size(A,1),1);u = u/norm(u);
end

if A ==0
    u = projectL0(u,param.k);v = projectL0(v,param.q);
elseif param.PSD
     u = conGradU02(A ,param.powerIter,param.k); v = u;
else
    pold=u'*A*v;
    for i = 1:param.powerIter
%        fprintf('Powerit %d %d..',length(find(u)),length(find(v)))
        v = A'*u;
%        fprintf('Powerit %d\n',length(find(v)))
        v = projectL0(v,param.q);
        v = v/norm(v);
        
        u = A*v;
        u = projectL0(u,param.k);
        u = u/norm(u);
        
        pnew=u'*A*v;
        err=abs((pold-pnew)/pold);
        if err<1e-8
            if debug
             fprintf(['nb iterations in truncatedPowerSeq ' num2str(i) '\n' ]);
            end
            break;
        end
        pold=pnew;
    end
    if debug==1 && i==param.powerIter
        if debug
        fprintf(['In truncatedPowerSeq exit error=' num2str(err) '\n']);
        end
    end
end


end