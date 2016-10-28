function [Sh,tracenorm,u,s] = ShrinkPSD(S,tau)
    S = .5*(S+S');
    [u,s] = eig(S);
    s=max(real(s),0);
    [ds,o] = sort(diag(s),'descend');
    u = u(:,o);
    s = s(o,o);
    r = sum(ds>tau);
    Sh = u(:,1:r)*diag(ds(1:r)-tau)*u(:,1:r)';
    Sh = real(Sh)+1e-16*eye(size(S,1)); %avoid numerical issues
%     eigs(Sh)
%     keyboard;
    tracenorm = sum(ds(1:r)-tau);
    u = u(:,1:r);
    s = s(1:r,1:r) - tau*eye(r);
    s = diag(s);
end