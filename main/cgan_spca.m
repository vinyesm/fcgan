function [Z ActiveSet hist param flag output] = cgan_spca(inputData,param,startingZ,ActiveSet)

% Active set algorithm to minimize f(Z) + lambda Omega_kq(Z), where f is a
% convex function with Lipschitz gradient and Omega_kq is the kq-trace norm
%  reference:
%% "Tight convex relaxations for sparse matrix factorization",
% by E. Richard, G. Obozinski, J.-P. Vert 2014.
% code: E. Richard, G. Obozinski, J.-P. Vert 2014.
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT
%
%% param contains
%
% param.f : the function to minimize, encoded with an integer param.f =
% 1,2,3 or 4:
%       (1) stands for L2 distance squared : f(Z) = .5 |Z-Y|^2
%       (2) stands for multitask :           f(Z) = .5 |XZ - Y|_F^2
%       (3) stands for quadratic regression  f(Z) = .5 | diag( X Z X' ) - y|_2^2
%       (4) stands for bilinear regression   f(Z) = .5 |X1 Z X2 - Y|_F^2
%
% param.k param.q : block size parameters
% param.lambda : lagrange multiplier
% param.cardfun : cardinality function. A vector of size size(Z,1)=size(Z,2)
%
%  ** OTHER PARAMETERS HAVE DEFAULT VALUES **
%
% param.PSD : true if the ouput is desired to be PSD, default value is false
% max_nb_main_loop : nb of runs of "Solve P_S", default = 100
% param.powerIter : nb runs power iteration, default = 100;
% param.stPtPowerIter : nb of different starting points for truncated power
% teration, default = 100
% param.epsStop : add a component to the working set when  varIJ >
% lambda*(1-param.epsStop), default = .1
%
% param.innerLoopIter : iterations of the prox-grad loop P_S default = 100
% param.niterPS : number of times an active component is picked and the
% relative subproblem solved, default = 200
%
% param.PSdualityEpsilon : tolerance for duality gap of subproblems,
% default : 1e-3

%% inputData contains
% inputData.X, inputData.Y X and Y used to compute the loss and the gradient
%
%% startingZ
% initial value of Z, for warm start, if not provided it is set to 0
%
%% ActiveSet
% use this only if initial Z is nonzero. It should be the ActiveSet object
% related to the starting Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT
% Z : the minimizer, a matrix of size m1 x m2
%
% ActiveSet : a list of elements containing
% ActiveSet.I : the sets suggested by the heuristic
% ActiveSet.J : the sets suggested by the heuristic (right)
% ActiveSet.U : the pseudo-singular vectors (of size k)
% ActiveSet.V : the pseudo-singular vectors (of size q)
% ActiveSet.Sigma : the pseudo-singular values
% ActiveSet.Z : the matrices of k x q such that Z = sum_i Z{i}
% ActiveSet.tracenorm : the tracenorms of the Z{i}'s (may be equal to 0)
% ActiveSet.fronorm : the frobenius square norms of the Z{i}'s (may be
% equal to 0), used to ad strong convexity to PS
% ActiveSet.k : Size of the Z{i}'s (kxk)
% ActiveSet.left ActiveSet.middle ActiveSet.right are build so that their product equals
% the solution:
% Z = ActiveSet.left*ActiveSet.middle*ActiveSet.right
%
%
% hist : a list of vectors containing the objective function values through
% iterations, the loss, penalty and duality gap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_store_in_asqp=0;

%%
if nargin < 3
    startingZ = set_default_Z(inputData,param);
    Z = startingZ;
    ActiveSet = {};
    ActiveSet.I = {};
    ActiveSet.U = {};
    ActiveSet.Sigma = {};
    ActiveSet.Z = {};
    ActiveSet.tracenorm = {};
    ActiveSet.fronorm = {};
    ActiveSet.k = {};    
    ActiveSet.atomsSupport = {};
    ActiveSet.alpha= [];
    ActiveSet.atom_count = 0;
    ActiveSet.max_atom_count_reached=0;
else
    Z = startingZ;
end

param = set_default_param(param);

obj  = [];
loss = [];
pen  = [];
dg   = [];
time =[];
dg_sup = [];
time_sup =[];
obj_sup = [];
nb_pivot=[];
active_var=[];
hist.varIJ=[];
hist.nzalphas=[];
hist.normalpha=[];
p=size(Z,1);
D=zeros(p,1);
output.time=0;

tic

max_nb_atoms = param.max_nb_atoms;      
max_nb_main_loop = param.max_nb_main_loop;
hist.time=toc;



%%
c = 1;
i = 0;



while c
    i = i+1;
    
    if ActiveSet.atom_count > max_nb_atoms
        break;
    end
    
    %% solve problem P_S
    
    if ~isempty(ActiveSet.I)
        if param.verbose==1
            fprintf('    solving PS..\n ')
        end
        if strcmp('asqp',param.opt)
            [Z, ActiveSet, hist_ps] = solve_ps_spca_asqp(Z, ActiveSet,param,inputData);
            if ~isempty(ActiveSet.alpha) && param.debug==1
                hist.nzalphas=[hist.nzalphas full(sum(ActiveSet.alpha>0))];
                hist.normalpha=[hist.normalpha full(sum(ActiveSet.alpha))];
                fprintf('   nz alphas=%d  |alpha|=%f  dual gap =%f\n',full(sum(ActiveSet.alpha>0)),full(sum(ActiveSet.alpha)), dg(end));
            end
            nb_pivot=[nb_pivot hist_ps.nb_pivot];
            active_var=[active_var hist_ps.active_var];
        end
        if strcmp('proxbcd',param.opt)
            [Z, D, ActiveSet, hist_ps] = solve_ps_spca_proxbcd(Z,D,ActiveSet,param,inputData);
        end
        
        obj = [obj hist_ps.obj];
        loss = [loss hist_ps.loss];
        pen = [pen hist_ps.pen];
        dg = [dg hist_ps.dg];
        time = [time hist_ps.time];
        dg_sup = [dg_sup hist_ps.dg_sup];
        time_sup = [time_sup hist_ps.time_sup];
        obj_sup = [obj_sup hist_ps.obj_sup];
    end
    
    %% get a new descent direction using truncated power iteration
    
    
    H = gradient(Z,inputData,param);

    
    if param.verbose==1
        fprintf('%d/%d   \n',i,max_nb_main_loop);
    end
    
    [u, kBest] = lmo_spca(-H,param);
    param.k=kBest;
    currI = find(u);
    
    
    
    %% verbose
    if param.verbose==1
        fprintf('   currI = ')
        for j=1:length(currI)
            fprintf('%d ',currI(j));
        end
        fprintf('\n');
        
        if(isempty(currI))
            fprintf('currI is empty\n');
        end
    end
    %%
    
    varIJ = norm(H(currI,currI));
    takenI= isInCell(currI,ActiveSet.I,cell2mat(ActiveSet.k)) ;
    hist.varIJ=[hist.varIJ varIJ];
    
    flag.var=varIJ;
    
    if param.verbose==1
        fprintf('   variance = %2.4e, thresh = %2.4e, length(currI)=%d\n',varIJ, param.lambda*(1+param.epsStop / kBest)* param.cardfun(kBest), length(currI))
    end
    
    
    if takenI
        fprintf('This support has already been added. Stopping\n\n');
        c=0;
    elseif varIJ > param.lambda*(1+param.epsStop / kBest)* param.cardfun(kBest)
        ActiveSet.I = [ActiveSet.I, currI];
        ActiveSet.U = [ActiveSet.U, u(currI)];
        ActiveSet.Sigma = [ActiveSet.Sigma, varIJ];
        ActiveSet.Z = [ActiveSet.Z, zeros(param.k,param.k)];
        ActiveSet.tracenorm = [ ActiveSet.tracenorm , 0];
        ActiveSet.k = [ActiveSet.k , kBest];
        ActiveSet.fronorm = [ ActiveSet.fronorm , 0];
    else
        c = 0;
    end
    c = i<max_nb_main_loop & c;
end

if param.debug==1
    if i>=max_nb_main_loop
        fprintf('\n max number of main loop iterations reached\n');
    end
end

% ActiveSet = postProcessFactors(ActiveSet,Z);
hist.obj = obj;
hist.loss = loss;
hist.pen = pen;
hist.dg = dg;
hist.time = time;
hist.dg_sup = dg_sup;
hist.time_sup= time_sup;
hist.obj_sup = obj_sup;
hist.nb_pivot= nb_pivot;
hist.active_var= active_var;

end

