function [ Z, ActiveSet, hist] = solve_ps_spca_asqp( Z,ActiveSet,param,inputData)
%Using Active Setet to solve (PS) problem


param_as.max_iter=1e3;
param_as.epsilon=1e-14;
param_as.debug_mode=false;
param_as.ws=true;

obj = zeros(1,param.niterPS);
pen =  zeros(1,param.niterPS);
loss =  zeros(1,param.niterPS);
dg =  zeros(1,param.niterPS);
time =  zeros(1,param.niterPS);
nb_pivot=zeros(1,param.niterPS);
active_var=zeros(1,param.niterPS);

p=size(Z,1);
cardVal=[];
lambdas=[];

if ActiveSet.atom_count>0
    if param.diag==1
        lambdas=param.lambda*[zeros(p,1);ones(ActiveSet.atom_count-p,1)];
        %precompute cardinality function values f(k)
        cardVal=sparse(p,1);
        for j=(p+1):ActiveSet.atom_count
            atom=ActiveSet.atoms(:,j);
            %     weight=param.cardfun(sum(atom~=0));
            weight=1; % modify when different weights for atoms
            cardVal=[cardVal;weight];
        end
    else
        lambdas=param.lambda*ones(ActiveSet.atom_count,1);
        for j=1:ActiveSet.atom_count
            atom=ActiveSet.atoms(:,j);
            %     weight=param.cardfun(sum(atom~=0));
            weight=1; % modify when different weights for atoms
            cardVal=[cardVal;weight];
        end
    end
    U=real(inputData.X1)*ActiveSet.atoms(:,1:ActiveSet.atom_count);
    G=U'*U;
    H=G.*G;
    % f= sum(U.*U)';
    f=diag(U'*inputData.Y*U);
    f=lambdas-f;
else
    U=[];
    G=[];
    H=[];
    f=[];
end

if param.debug
    alphaSparsity=[];
    diffalphas=[];
end


new_atom_added=false;
first_pass=1;
i=1;
count=1;
cont=true;

while cont
    %% get a new atom ui  [si=ui'*ui] in ActiveSet.atoms  from some C_I I
    % in ActiveSet.I to add to the collection
    % current point is Z
    
    if ActiveSet.atom_count>param.max_nb_atoms
        display('maximum number of atoms added');
        break;
    end
    
    if ~first_pass
        
        if new_atom_added
            alpha0=[ActiveSet.alpha;0];
        else
            alpha0=ActiveSet.alpha;
        end
        
        % active-set
        [alph,Jset,npiv]=asqp(H+1e-14*eye(ActiveSet.atom_count),-f,alpha0,param_as,1);
        
        new_atom_count=sum(Jset);
        ActiveSet.alpha=alph(Jset);
        ActiveSet.atom_count=new_atom_count;
        ActiveSet.atoms(:,1:new_atom_count)=ActiveSet.atoms(:,Jset);
        U=U(:,Jset);
        f=f(Jset);
        cardVal=cardVal(Jset);
        G=U'*U;
        H=G.*G;
        
        %% Update ActiveSet and Z
        Z=zeros(p);
        nz=find(ActiveSet.alpha>1e-15);
        for j=nz'
            u=ActiveSet.atoms(:,j);
            Z=Z+ActiveSet.alpha(j)*(u*u');
        end
        
        %% Compute objective, loss, penalty and duality gap
        if param.sloppy==0 || (param.sloppy~=0 && mod(count,10)==1)
            [loss(i),pen(i),obj(i),dg(i),time(i)]=get_val_spca_asqp(ActiveSet,inputData,param,cardVal);
            nb_pivot(i)=npiv;
            active_var(i)= sum(ActiveSet.alpha>0);
            cont = (dg(i)>param.PSdualityEpsilon) && count< param.niterPS;
            i=i+1;
        end
    end
    
    %% get new atom
    if cont
        [new_i, new_val, maxval]=get_new_atom_spca(Z,ActiveSet,param,inputData);
        ActiveSet.atom_count = ActiveSet.atom_count +1;
        ActiveSet.max_atom_count_reached=max(ActiveSet.max_atom_count_reached,ActiveSet.atom_count);
        ActiveSet.atomsSupport=[ActiveSet.atomsSupport new_i];
        anew=sparse(new_i,ones(length(new_i),1),new_val,p,1);
        ActiveSet.atoms(:,ActiveSet.atom_count)=anew;
        
        unew=real(inputData.X1)*anew;
        if isempty(U)
            vnew=[];
        else
            vnew=U'*unew;
            vnew=vnew.*vnew;
        end
        U=[U unew];
%         G2=U'*U;        
%         Gold=G;
%         G=zeros(ActiveSet.atom_count);
%         G(1:ActiveSet.atom_count-1,1:ActiveSet.atom_count-1)=Gold;
%         G(1:ActiveSet.atom_count-1,end)=vnew;
%         G(end,1:ActiveSet.atom_count-1)=vnew';
%         G(end,end)=unew'*unew;
        Hold=H;
        H=zeros(ActiveSet.atom_count);
        H(1:ActiveSet.atom_count-1,1:ActiveSet.atom_count-1)=Hold;
        H(1:ActiveSet.atom_count-1,end)=vnew;
        H(end,1:ActiveSet.atom_count-1)=vnew';
        H(end,end)=(unew'*unew)^2;
%         H2=G2.*G2;
        weight=1;
        %     f= [f; param.lambda*weight-norm(unew)^2];
        f= [f; param.lambda*weight-unew'*inputData.Y*unew];
        cardVal=[cardVal;weight];
        
        new_atom_added=true;
        first_pass=false;
        
        if maxval<0
            error('\nNegative directional derivative d=%f\n',maxval);
        end
    end
    
    count=count+1;
end

%redimensioning arrays
i=i-1;
hist.obj = obj(1:i);
hist.pen =  pen(1:i);
hist.loss =  loss(1:i);
hist.dg =  dg(1:i);
hist.time = time(1:i);
hist.obj_sup=obj(1);
hist.dg_sup=dg(1);
hist.time_sup=time(1);
hist.nb_pivot=nb_pivot(1:i);
hist.active_var=active_var(1:i);

if count>param.niterPS
    fprintf('maximum number of Ps iteration reached, duality gap=%f\n',dg(end));
end


if param.debug
    figure(10);clf;
    subplot(1,4,1);
    plot(obj);
    title('objective');
    subplot(1,4,2);
    semilogy(dg);
    title('duality gap');
    subplot(1,4,3)
    plot(alphaSparsity);
    title('alpha nz');
    subplot(1,4,4)
    plot(diffalphas);
    title('diff alphas');
    keyboard;
end

end


