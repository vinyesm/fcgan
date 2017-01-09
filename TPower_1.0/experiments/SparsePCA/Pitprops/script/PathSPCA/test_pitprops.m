% Test code on various examples
clear all;

load ..\..\data\pitprops.mat
S = corr_matrix_pitprops;

%% -------------------------------- PathSPCA for Sparse PCA -------------------------------
cardinality_vec = [7 2 1 1 1 1];

disp('PathSPCA: Computing Sparse Loadings with Cardinality:');
disp(cardinality_vec);

tic;
[V]= PartialPathCov_Multiple(S, cardinality_vec);%FullPathCov(S);
t_pathspca = toc;

variance_rate = vars_adj(V, S);

fprintf('PathSPCA: Propotion of explained variance: %f...\n', variance_rate);
% bnds_rate = bnds / trace(covcolon);

% save results
PathSPCA_Result = [];
PathSPCA_Result.CPUTime = t_pathspca;
PathSPCA_Result.variance_rate = variance_rate;
PathSPCA_Result.loadings = V;

save ..\..\result\PathSPCA_Result.mat 'PathSPCA_Result';