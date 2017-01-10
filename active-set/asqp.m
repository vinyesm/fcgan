function [c,A,nbpivot]=asqp(Q,b,c0,param,new_atom_added)
%
% SOLVING QUADRATIC PROBLEM WITH ACTIVE SET, WITH WARM START
% min_c 0.5*c'Qc-c'*b s.t. c \geq 0

% INPUTS:
% Q,b,c0           = quadratic problem parameters and intitial point c0
% new_atom_added   = weather or not a new atom has been added to the active set
% param.max_iter   = max nb of iterations
% param.epsilon    = tolerance parameter
% param.debug_mode = true/false
%
% OPTIONAL INPUTS
% param.max_nb_iter    = default 500 = maximum number of iterations allowed default = 500;
% param.max_nb_atoms   = default 500 = maximum number of atoms (used to allocate memory)
% param.epsStop        = default 1e-5 = tolerance parameter
% param.flag_no_design = when no design matrix more optimized storage.
% param.debug          = false = to debug;
% param.debug_asqp     = false = to debug inside asqp
%
% OUTPUTS :
% c       = final result
% A       = vector of 0/1 indicating active atoms
% nbpivot = final nb of pivots, i.e. full_steps+drop_steps
%
%%%%%%%%%%%%
%
% Active-set algorithm from
% Numerical Optimization - J.Nocedal & S.J.Wright
%
% Marina Vinyes and Guillaume Obozinski, 2016
% %%%%%%%%%%%%


max_iter=param.max_iter;
epsilon=param.epsilon;
debug_mode=param.debug_mode;
debug_mode=1;

tol=1e-14;
t=size(c0,1);
c=full(c0);
g=Q*c-b;
A=(c>tol); % active atoms
if param.ws
    if new_atom_added,
        A(end)=true; % adding the last one
        if g(end)>0,
            error('The new atom direction is not a descent direction')
        end
    end
else
    if(any(g<-tol))
        [~,j]=min(g);
        A(j)=true;
    else
        error('No direction added')
    end
end

hist.norm_g=zeros(t,max_iter);
nb_drop_steps=0;
nb_full_steps=0;

if debug_mode,
    hist.obj=zeros(1,max_iter);
    hist.J=zeros(t,max_iter);
    hist.c=zeros(t,max_iter);
    hist.r=zeros(1,max_iter);
    obj_old=0.5*c'*Q*c-c'*b;
end

iter=1;
while(iter<=max_iter)
    %% Compute new candidate solution
    J=A|(g<0);
    if norm(g(J))<epsilon,
        break;
    else
        hist.norm_g(iter)=norm(g(A));
    end
    d=zeros(t,1);
    if debug_mode && (rcond(Q(A,A))<1e-12),
        display('asqp : Warning, Hessian badly conditioned\n');
        %keyboard;
    end
    d(A)=Q(A,A)\b(A);
    %% Progress until active set reduces
    if (~all(d(A)>=0)), % Drop step
        K=A & (c-d>0);
        [tau,i_remove]=min(c(K)./(c(K)-d(K)));
        if (tau>1),
            error('problem with tau');
        end
        I_K=find(K);
        A(I_K(i_remove))=false;
        c=c+tau*(d-c);
        c(I_K(i_remove))=0;
        c(c<tol)=0;
        g=Q*c-b;
        A=A & (c>0);
        nb_drop_steps=nb_drop_steps+1;
        %fprintf('.');
    else % Full step
        c=d;
        g=Q*c-b;
        nb_full_steps=nb_full_steps+1;
        %fprintf('+')
        if param.ws && nb_full_steps>10,
            break;
        end
        if(any(g<-tol & ~A))
            [~,j]=min(g.*(~A));
            A(j)=true;
        end
    end
    if debug_mode,
        hist.obj(iter)=0.5*c'*Q*c-c'*b;
        if hist.obj(iter)>obj_old+tol,
            error('obj increases in asqp');
        end
    end
    %% Test to increase active set
    if any(A & c==0 & g>=0)
        error('wrong feature');
    end
    if debug_mode,
        obj_old=hist.obj(iter);
        hist.J(:,iter)=A;
        hist.c(:,iter)=c;
        if(all(d(A)>=0)),
            hist.r(iter)=0;
        else
            hist.r(iter)=I_K(i_remove);
        end
    end
    iter=iter+1;
end

nbpivot=nb_full_steps+nb_drop_steps;
%fprintf('\n');

if debug_mode,
    hist.obj=hist.obj(1:min(iter,max_iter)-1);
    hist.norm_g=hist.norm_g(1:min(iter,max_iter));
    if any(diff(hist.obj)>tol),
        error('objective increases in asqp');
    end
    if iter>max_iter,
        figure(15);
        plot(hist.obj);
        keyboard;
    end
end



