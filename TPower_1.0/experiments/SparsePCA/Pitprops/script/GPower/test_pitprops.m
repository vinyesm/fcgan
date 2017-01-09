clear all; 
load ..\..\data\pitprops.mat

S = corr_matrix_pitprops;
A = chol(S);
% Set up Objective Function

%% -------------------------------- GPower for Sparse PCA -------------------------------

m=6;                                % Number of components:

% ***** Single-unit algorithms *****
gamma=0.4*ones(1,m);                % sparsity weight factors -one for each component - 
   
tic;
Z1=GPower(A,gamma,m,'l1',0);        
t_cpu = toc;
variance_rate = vars_adj(Z1,S);

fprintf('GPower: Propotion of explained variance: %f...\n', variance_rate);

% save results
GPower_Result = [];
GPower_Result.CPUTime = t_cpu;
GPower_Result.variance_rate = variance_rate;
GPower_Result.loadings = Z1;

save ..\..\result\GPower_Result.mat 'GPower_Result';