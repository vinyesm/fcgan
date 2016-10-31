function [x, as, hist, iter] = cgan_lgl_constrained(D,y,param)
%
% usage [x, as, hist, iter] = cgan_lgl(D,y,param)
%
% FUNCTION PERFORMS A MINIMIZATION OF THE FORM:
% min_x 0.5|| y - D x ||^2  s.t.  ||x||_A \leq rho
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
param_as.max_iter= 100;
param_as.epsilon= 1e-10;
param_as.rho=param.rho;

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
rho=param.rho;

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
                    b=as.atoms(:,1:atom_count)'*y;
                else
                    b=as.DA(:,1:atom_count)'*y;
                end
                if atom_count==1,
                    coeffs=rho;
                    Jset=1;
                    npiv=0;
                else
                    %                     [alpha_qpas,Jset,npiv]=asqp_new_constrained2(H(1:atom_count,1:atom_count)+1e-12*eye(atom_count),b,alpha0,100,1e-10,debug_asqp,new_atom_added,rho2);
                    [coeffs,Jset,npiv]=asqp_constrained(H(1:atom_count,1:atom_count)+1e-12*eye(atom_count),b,coeffs_ws,param_as,new_atom_added);
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
                gamma=2/(iter+10);
                if new_atom_added,
                    coeffs=[(1-gamma)*as.coeffs;gamma*rho];
                else
                    coeffs=(1-gamma)*as.coeffs;
                end
            case 'FWls'
                if new_atom_added,
                    if param.flag_no_design,
                        fw_dir=rho*as.atoms(:,atom_count)-x;
                        gamma=min(max((-g'*fw_dir)./(sum(fw_dir.^2)),0),1);
                    else
                        fw_dir=rho*as.DA(:,atom_count)-Dx;
                        gamma=min(max((-(Dx-y)'*fw_dir)./(sum(fw_dir.^2)),0),1);
                    end
                    coeffs=[(1-gamma)*as.coeffs;gamma*rho];
                else
                    error('no new atom\n');
                end
            case 'FWpw'
                if new_atom_added,
                    if iter==1
                        coeffs=rho;
                    else
                        if param.flag_no_design,
                            score=g'*as.atoms(:,1:atom_count);
                            score(as.coeffs<1e-12)=-Inf;
                            [away_val,away_i]=max(score);
                            fw_dir=rho*(as.atoms(:,atom_count)-as.atoms(:,away_i));
                            maxgamma=as.coeffs(away_i)/rho2;
                            gamma=min(max((-((Dx-y)'*fw_dir))./(sum(fw_dir.^2)),0),maxgamma);
                        else
                            score=g'*as.atoms(:,1:atom_count);
                            score(as.coeffs<1e-12)=-Inf;
                            [away_val,away_i]=max(score);
                            fw_dir=rho*(as.DA(:,atom_count)-as.DA(:,away_i));
                            maxgamma=as.coeffs(away_i)/rho;
                            gamma=min(max((-((Dx-y)'*fw_dir))./(sum(fw_dir.^2)),0),maxgamma);
                        end
                        coeffs=[as.coeffs;gamma*rho];
                        coeffs(away_i)=coeffs(away_i)-gamma*rho;
                        if coeffs(away_i)<-1e-10,
                            error('coeffs(away_i) negative');
                        elseif coeffs(away_i)<1e-12,
                            coeffs(away_i)=0;
                        end
                    end
                else
                    error('no new atom\n');
                end
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
    
    pen(iter)=tau;
    obj(iter)=loss(iter);
    
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
    maxval=rho*maxval;
    
    if maxval>0,
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
    if param.flag_no_design,
        dg(iter)=g'*x+maxval;
    else
        dg(iter)=g'*x+maxval;
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





