function [x, as, hist, iter] = cgan_lgl(D,y,param)
%
% usage [x, as, hist, iter] = cgan_lgl(D,y,param)
%
% FUNCTION PERFORMS A MINIMIZATION OF THE FORM:
% min_x 0.5|| y - D x ||^2  + lambda * ||x||_A 
%
%
% INPUTS:
% y, D         = observarion and data matrix
% param.lmo    = function handle to pick an atom. one of its inputs must be the gradient.
% param.lambda = regularization parameter for the atomic norm
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
% x           = final result
% as.atoms    = final atomic set
% as.DA       = D*atoms
% as.coeffs   = final coefficients
% hist.obj    = objective function value
% hist.dg     = duality gap
% hist.loss   = loss  value
% hist.pen    = penality value   
% hist.tt     = time taken since start of the iterations
% iter        = nb of iterations in total (for column generation, the nb of calls to active-set)
%
%%%%%%%%%%%%
% Marina Vinyes and Guillaume Obozinski, 2016
% %%%%%%%%%%%%


param_as.debug_mode= false;

param_as.epsilon= 1e-10;
param_as.ws=param.ws;
if param.ws
    param_as.max_iter= 1000;
else
    param_as.max_iter= 100;
end

if param.flag_no_design,
    p=length(y);
    n=p;
    if ~isempty(D),
        error('No design is expected. Set D=[]');
    end
    x=zeros(p,1);
else
    [n,p]=size(D);
    x=zeros(p,1);
    Dx=zeros(n,1);
end

max_nb_atoms=param.max_nb_atoms;
max_nb_iter=param.max_nb_iter;
lambda=param.lambda;
rho=0.5*sum(y.^2)/lambda;

as.atoms=sparse([],[],[],p,max_nb_atoms,max_nb_atoms*param.K);
as.DA=sparse([],[],[],n,max_nb_atoms,max_nb_atoms*param.K);
H=sparse([],[],[],max_nb_atoms,max_nb_atoms,max_nb_atoms*max_nb_atoms);
as.coeffs=[];

if param.debug
    maxvals=[];
end

dg=zeros(max_nb_iter,1);
obj=zeros(max_nb_iter,1);
loss=zeros(max_nb_iter,1);
pen=zeros(max_nb_iter,1);
tt=zeros(max_nb_iter,1);
nb_pivot=zeros(max_nb_iter,1);
active_var=zeros(max_nb_iter,1);

npiv=0;
max_atom_count_reached=0;
iter=0;
atom_count=0;

tic

while(iter<max_nb_iter),
    if iter>0,
        if atom_count==0,
            error('all atoms have been thrown away');
        end
        if new_atom_added,
            coeffs_ws=[as.coeffs;0];
        else
            coeffs_ws=as.coeffs;
        end
        switch param.method,
            case 'asqp'
                if param.flag_no_design,
                    H=as.atoms(:,1:atom_count)'*as.atoms(:,1:atom_count);
                    b=as.atoms(:,1:atom_count)'*y-lambda*ones(atom_count,1);
                else
                    b=as.DA(:,1:atom_count)'*y-lambda*ones(atom_count,1);
                end
                if param.ws
                    [coeffs,Jset,npiv]=asqp(H(1:atom_count,1:atom_count)+1e-10*eye(atom_count),b,coeffs_ws,param_as,new_atom_added);
                else
                    [coeffs,Jset,npiv]=asqp(H(1:atom_count,1:atom_count)+1e-10*eye(atom_count),b,zeros(atom_count,1),param_as,new_atom_added);
                end
                % Hard threshold small negative values
                smallValues=find(coeffs<0); % for numerical issues
                coeffs(smallValues)=zeros(length(smallValues),1);
                %manage H and the set of atoms
                atom_count=sum(Jset);
                H=H(Jset,Jset);
                as.atoms(:,1:atom_count)=as.atoms(:,Jset);
                as.atoms(:,(atom_count+1):end)=0;
                as.DA(:,1:atom_count)=as.DA(:,Jset);
                as.DA(:,(atom_count+1):end)=0;
                coeffs=coeffs(Jset,1);               
            case 'FW'
                gamma=2/(iter+1000);
                if new_atom_added,
                    coeffs=[(1-gamma)*as.coeffs;gamma*rho];
                else
                    coeffs=(1-gamma)*as.coeffs;
                end
            case 'FWls'
                if new_atom_added,
                    if param.flag_no_design,
                        fw_dir=rho*as.atoms(:,atom_count)-x;
                        gamma=min(max((-(g'*fw_dir+lambda*rho-pen(iter)))./(sum(fw_dir.^2)),0),1);
                    else
                        fw_dir=rho*as.DA(:,atom_count)-Dx;
                        gamma=min(max((-((Dx-y)'*fw_dir+lambda*rho-pen(iter)))./(sum(fw_dir.^2)),0),1);
                        
                    end
                    coeffs=[(1-gamma)*as.coeffs;gamma*rho];
                else
                    if param.flag_no_design,
                        fw_dir=-x;
                        gamma=min(max((-(g'*fw_dir-pen(iter)))./(sum(fw_dir.^2)),0),1);
                    else
                        fw_dir=-Dx;
                        gamma=min(max((-((Dx-y)'*fw_dir-pen(iter)))./(sum(fw_dir.^2)),0),1);
                    end
                    coeffs=(1-gamma)*as.coeffs;
                end
            case 'FWpw'
                if new_atom_added,
                    if param.flag_no_design,
                        score=g'*as.atoms(:,1:atom_count);
                        score(as.coeffs<1e-12)=-Inf;
                        [away_val,away_i]=max(score);
                        if away_val>=-lambda-1e-12,
                            fw_dir=rho*(as.atoms(:,atom_count)-as.atoms(:,away_i));
                            gamma=min(max((-(g'*fw_dir))./(sum(fw_dir.^2)),0),as.coeffs(away_i)/rho);
                        else
                            away_i=0;
                            fw_dir=rho*as.atoms(:,atom_count);
                            gamma=min(max(-(g'*fw_dir+lambda*rho)./(sum(fw_dir.^2)),0),1-tau/rho);
                        end
                    else
                        score=g'*as.atoms(:,1:atom_count);
                        score(as.coeffs<1e-12)=-Inf;
                        [away_val,away_i]=max(score);
                        if away_val>=-lambda-1e-12,
                            fw_dir=rho*(as.DA(:,atom_count)-as.DA(:,away_i));
                            gamma=min(max((-((Dx-y)'*fw_dir))./(sum(fw_dir.^2)),0),as.coeffs(away_i)/rho);
                        else
                            away_i=0;
                            fw_dir=rho*(as.DA(:,atom_count));
                            gamma=min(max(-((Dx-y)'*fw_dir+lambda*rho)./(sum(fw_dir.^2)),0),1-tau/rho);
                        end
                    end
                    if away_val>=-lambda-1e-12,
                        coeffs=[as.coeffs;gamma*rho];
                        coeffs(away_i)=coeffs(away_i)-gamma*rho;
                        if coeffs(away_i)<-1e-10,
                            error('coeffs(away_i) negative');
                        elseif coeffs(away_i)<1e-12,
                            coeffs(away_i)=0;
                        end
                    else
                        coeffs=[as.coeffs;gamma*rho];
                    end
                else % Forward step is towards the origin
                    if param.flag_no_design,
                        score=g'*as.atoms(:,1:atom_count);
                        score(as.coeffs<1e-12)=-Inf;
                        [away_val,away_i]=max(score);
                        if away_val>=-lambda-1e-12,
                            fw_dir=-rho*as.atoms(:,away_i);
                            gamma=min(max((-(g'*fw_dir-lambda*rho))./(sum(fw_dir.^2)),0),as.coeffs(away_i)/rho);
                        else
                            away_i=0;
                            display('Away step and forward steps are both towards the origin');
                            break;
                        end
                    else
                        score=g'*as.atoms(:,1:atom_count);
                        score(as.coeffs<1e-12)=-Inf;
                        [away_val,away_i]=max(score);
                        if away_val>=-lambda-1e-12,
                            fw_dir=-rho*as.DA(:,away_i);
                            gamma=min(max((-((Dx-y)'*fw_dir-lambda*rho))./(sum(fw_dir.^2)),0),as.coeffs(away_i)/rho);
                        else
                            away_i=0;
                            display('Away step and forward steps are both towards the origin');
                            break;
                        end
                    end
                    if away_val>=-lambda-1e-12,
                        coeffs=as.coeffs;
                        coeffs(away_i)=coeffs(away_i)-gamma*rho;
                        if coeffs(away_i)<-1e-10,
                            error('coeffs(away_i) negative');
                        elseif coeffs(away_i)<1e-12,
                            coeffs(away_i)=0;
                        end
                    else
                        display('Away step and forward steps are both towards the origin');
                        break;
                    end
                end
                if gamma==0,
                    error('Stuck!');
                    break;
                end
%             case 'FWfully'
%                 coeffs = solve_fully_cvx(y,as.DA(:,1:atom_count),rho,as.coeffs);
            otherwise
                error('Unknown method');
        end    
        
        %% Compute the current solution
        as.coeffs=sparse(coeffs);
        x=as.atoms(:,1:atom_count)*as.coeffs;
        if ~param.flag_no_design,
            Dx=as.DA(:,1:atom_count)*as.coeffs;
            x=as.atoms(:,1:atom_count)*as.coeffs;
        end        
    end
    
    iter=iter+1;
    
    %% Compute gradient
    if param.flag_no_design,
        g=x-y;
    else
        g=D'*(Dx-y);
    end  
    
    
    %% Compute objective, loss and penalty
    
    tau=sum(as.coeffs);
    if param.flag_no_design,
        loss(iter)=0.5*sum(g.^2);
    else
        loss(iter)=0.5*sum((Dx-y).^2);
    end
    
    pen(iter)=lambda*tau;
    obj(iter)=loss(iter)+pen(iter);

    if iter>1,
        if obj(iter)>obj(iter-1)
        end
    end
    
    %%  Store nb of pivot and nb of active variables
    
    if strcmp(param.method,'asqp')
        nb_pivot(iter)=npiv;
        active_var(iter)=sum(as.coeffs>0);
    end
    
    %% Get new atom    
    
    [maxval,new_atom]=feval(param.lmo,-g,param);
    
    if maxval>lambda,
        A=1:atom_count;
        atom_count=atom_count+1;
        max_atom_count_reached=max(max_atom_count_reached,atom_count);
        % Inserting the new atom
        if param.flag_no_design,
            as.atoms(:,atom_count)=new_atom;
        else
            as.atoms(:,atom_count)=new_atom;
            as.DA(:,atom_count)=D*new_atom;
            % rank one update of H
            if strcmp(param.method,'asqp')
                vA=as.DA(:,A)'*as.DA(:,atom_count);
                aA= as.DA(:,atom_count)'*as.DA(:,atom_count);
                %             HA=H(A,A);
                %             HAvA=HA*vA;
                %             beta=(aA-vA'*HA*vA)^-1;
                %             H(A,A)=HA+beta*(HAvA*HAvA');
                %             H(atom_count,atom_count)=beta;
                %             H(A,atom_count)=-beta*HAvA;
                %             H(atom_count,A)=-beta*HAvA';
                H(A,atom_count)=vA;
                H(atom_count,A)=vA';
                H(atom_count,atom_count)=aA;
            end
        end
        if full(new_atom)'*g>0,
            error('new atom wrong');
        end
        new_atom_added=true;
    else
        new_atom_added=false;
        switch param.method,
            case 'asqp',
                error('In ASQP an atom should be added at each step. Debug needed');
            otherwise
                %disp('Forward step towards origin');
        end
    end
    
    if atom_count>max_nb_atoms,
        error('max number of atoms reached. Either atom storage is too small or atom dropping is not working correctly');
    end
    
    
    %% Compute duality gap
    
    c=min(1,lambda./maxval);
    if param.flag_no_design,
        dg(iter)=0.5*(1-c)^2*sum(g.^2)+lambda*tau+c*g'*x;
    else
        dg(iter)=0.5*(1-c)^2*sum((Dx-y).^2)+lambda*tau+c*(Dx-y)'*Dx;
    end
    
    tt(iter)=toc;
    
    if dg(iter) <= param.epsStop,
        if param.debug,
            disp('Terminating successfully with small duality gap');
        end
        break;
    end
    
    %% Debug mode
    
    if param.debug && iter>1,
        maxvals=[maxvals maxval];
    end
    
end

if iter>= max_nb_iter,
    disp('Max number of iterations reached in solve_fw');
end


%% Hist outputs

iter=min(iter,max_nb_iter);
burn_in=min(1,iter);

hist.obj = obj(burn_in:iter);
hist.pen =  pen(burn_in:iter);
hist.loss =  loss(burn_in:iter);
hist.dg =  dg(burn_in:iter);
hist.tt =  tt(burn_in:iter);
if strcmp(param.method,'asqp')
    hist.nb_pivot=nb_pivot(burn_in:iter);
    hist.active_var= active_var(burn_in:iter);
end

%% Debug figures
if param.debug
    
    figure(10);clf;
    subplot(2,3,1);
    plot(burn_in:iter,hist.obj,'b.');
    title('objective');
    subplot(2,3,2);
    plot(burn_in:iter,hist.loss,'b.');
    title('loss');
    subplot(2,3,3);
    plot(burn_in:iter,hist.pen,'b.');
    title('\lambda.NormUpBd');
    subplot(2,3,4);
    semilogy(burn_in:iter,hist.dg,'b-');
    title('duality gap');
    subplot(2,3,5)
    plot(maxvals);
    title('maxvals');  
    
    
    figure(11)
    if param.flag_no_design,
        imagesc(as.atoms(:,1:atom_count));
        xlabel('atom index');
        title('Atoms selected');
    else
        imagesc(as.atoms(:,1:atom_count));
        xlabel('atom index');
        title('D.Atoms selected');
    end
    if strcmp(param.method,'asqp')
        total_piv=cumsum(hist.nb_pivot);
        total_piv=total_piv(burn_in:iter);
        figure(15);clf;
        subplot(2,2,1);
        bar(hist.nb_pivot);
        xlim([0 length(total_piv)])
        pbaspect([1,1,1])
        ylabel('#pivots');
        xlabel('iterations');
        title('#pivots per Active-Set call');
        subplot(2,2,2);
        bar(total_piv);
        xlim([0 length(total_piv)])
        pbaspect([1,1,1])
        xlabel('iterations');
        ylabel('#total pivots');
        title('total #pivots');
        subplot(2,2,3);
        bar(hist.active_var,'FaceColor',[.6 0 0]);
        pbaspect([2,1,1])
        xlabel('iterations');
        ylabel('#active variables');
        %legend('active', 'non-active');
        title('#active variables during iterations');
    end
    keyboard;
    
end





