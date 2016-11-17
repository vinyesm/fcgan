function [c,x_ws,gamma_ws,A,nbpivot]=as_tr_l1(y,atoms,c0,x0,gamma0,param)

% [x,d,gamma] solve subproblem with bcmm
% s2 = gamma
% test <s2,ai><= mu forall i in Jc
% J <= J U {imax}, adding most violated constraint


debug_mode=1;

max_iter=param.max_iter;
% epsilon=param.epsilon;
lambda=param.lambda;
mu=param.mu;


param_bcmm.debug_mode=1;
param_bcmm.max_iter=10;
% param_bcmm.rho=1;
param_bcmm.lambda=lambda;
param_bcmm.mu=mu;

tol=1e-14;
t=size(c0,1);
c=full(c0);
x=x0;
A=(c>tol); % active atoms
A(end)=true;


hist.norm_g=zeros(t,max_iter);
nb_drop_steps=0;
nb_full_steps=0;

if debug_mode,
    hist.obj=zeros(1,max_iter);
    hist.J=zeros(t,max_iter);
    hist.c=zeros(t,max_iter);
    hist.r=zeros(1,max_iter);
    obj_old=0.5*norm(y-atoms*c,'fro')^2+lambda*sum(abs(x))+mu*sum(c);
end

gamma_ws=y;
x_ws=zeros(size(y,1),1);
d_ws=zeros(length(c0),1);
d=d_ws;

iter=1;
fprintf('as\n');

while(iter<=max_iter)
    %% Compute new candidate solution
%     J=A|(g<0);
%     if norm(g(J))<epsilon,
%         break;
%     else
%         hist.norm_g(iter)=norm(g(A));
%     end
%     d=zeros(t,1);
    d_ws=d;
    [d(A), x_ws, gamma_ws, r]=bcmm_tr_l1(y,atoms(:,A),d_ws(A),x_ws, gamma_ws, param_bcmm);
    %% Progress until active set reduces
    keyboard;
    K=find(~A);
    res=0;
    
    %test
    A0=A&(d<1e-12);
    Aind=find(A0);
    if ~isempty(K)
    dotprods= sum(bsxfun(@times,atoms(:,A0),gamma_ws)); 
    [res, ind]=max(dotprods);
    if res>mu
        error('problem with conditions');
    end
    end
    
    if ~isempty(K)
    dotprods= sum(bsxfun(@times,atoms(:,~A),gamma_ws));    
    [res, i_remove]=max(dotprods);
    end
        
%     keyboard;
    if ~isempty(K) && (res>mu), % Drop step
        A(K(i_remove))=false;
        nb_drop_steps=nb_drop_steps+1;
        %fprintf('.');
    else % Full step
        c=d;
        x=x_ws;
        nb_full_steps=nb_full_steps+1;
        break;
        %fprintf('+')
        if param.ws && nb_full_steps>10,
            break;
        end
%         if(any(g<-tol & ~A))
%             [~,j]=min(g.*(~A));
%             A(j)=true;
%         end
    end
    if debug_mode,
        hist.obj(iter)=0.5*norm(y-atoms*c,'fro')^2+lambda*sum(abs(x))+mu*sum(c);
%         if hist.obj(iter)>obj_old+tol,
%             error('obj increases in asqp');
%         end
    end
    %% Test to increase active set
%     if any(A & c==0 & g>=0)
%         error('wrong feature');
%     end
    if debug_mode,
        obj_old=hist.obj(iter);
    end
    iter=iter+1;
end

nbpivot=nb_full_steps+nb_drop_steps;
%fprintf('\n');

if debug_mode,
    hist.obj=hist.obj(1:min(iter,max_iter)-1);
    if any(diff(hist.obj)>tol),
        error('objective increases in asqp');
    end
    if iter>max_iter,
        fprintf('max number of iterations in as_tr_l1\n');
        figure(15);
        plot(hist.obj);
%         keyboard;
    end
end



