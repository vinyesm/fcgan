function [d, x, gamma,Jset, r]=bcmm_tr_l1(y,A,d0,x0, gamma0, param)

% solves min |y-Ad|^2+lambda * |x|_1 + mu* 1'*d s.t. Ad=x and d>=0

debug_mode=1;
max_iter=50;

x=x0;
d=d0;
gamma=gamma0;
lambda=param.lambda;
mu=param.mu;
rho=0.1*lambda;
H=(A'*A);
b0=A'*y-mu;
% keyboard;

param_as.debug_mode= false;
param_as.epsilon= 1e-10;
param_as.ws=1;
param_as.max_iter= 100;

if debug_mode
    obj=zeros(1,max_iter);
    lag=zeros(1,max_iter);
    eq=zeros(1,max_iter);
    dg=zeros(1,max_iter);
end

r=0;
alpha=rho;
while r<max_iter
    
    %update on gamma,x
    r=r+1;
    %alpha=rho/sqrt(100+r);    
    gamma=gamma+alpha*(x-A*d);
    x=soft_threshold(A*d-gamma/rho,lambda/rho);
%     keyboard;
    
    if debug_mode
        obj(r)=0.5*norm(y-A*d,'fro')^2+rho/2*norm(x-A*d,'fro')^2+lambda*sum(abs(x))+mu*sum(d);
        lag(r)=0.5*norm(y-A*d,'fro')^2+rho/2*norm(x-A*d,'fro')^2+lambda*sum(abs(x))+mu*sum(d) + dot(gamma,x-A*d);
        eq(r)=norm(x-A*d,'fro')^2;
        dg(r)=-dot(gamma,x-A*d);
%         pres(r)=
%         dres(r)=
    end
        
    %update on gamma,d
    r=r+1;
    %alpha=rho/sqrt(100+r);    
    %gamma=gamma+alpha*(x-A*d);
    b=A'*(rho*x+gamma+y)-mu;
    [d,Jset,npiv]=asqp(H+1e-10*eye(size(H,1)),b,d,param_as,0);
    
    if debug_mode
        obj(r)=0.5*norm(y-A*d,'fro')^2+rho/2*norm(x-A*d,'fro')^2+lambda*sum(abs(x))+mu*sum(d);
        lag(r)=0.5*norm(y-A*d,'fro')^2+rho/2*norm(x-A*d,'fro')^2+lambda*sum(abs(x))+mu*sum(d) + dot(gamma,x-A*d);
        eq(r)=norm(x-A*d,'fro')^2;
        dg(r)=-dot(gamma,x-A*d);
    end
end

x=A*d;
Jset=(abs(d)>1e-15);

if debug_mode
    figure(10);clf;
    subplot(1,4,1);
    plot(obj,'.');
    title('augmented obj');
    subplot(1,4,2);
    plot(lag,'.');
    title('augmented lag');
    subplot(1,4,3);
    plot(eq,'.');
    title('equality constraint')
    subplot(1,4,4);
    semilogy(dg,'.');
    title('dg ')
    keyboard
end


end