clear all; 
load ..\..\data\pitprops.mat

A = corr_matrix_pitprops;

%% -------------------------------- TPower for Sparse PCA -------------------------------
% parameters setting
cardinality_vec = [7 2 1 1 1 1];
options = [];
options.cardinality_vec = cardinality_vec;
options.OptTol = 1e-8;
options.MaxIter = 500;
options.verbose = 0;
options.initType = 2;

disp('TPower: Computing Sparse Loadings with Cardinality:');
disp(cardinality_vec);

% run TPower 
tic;
U = TPower_SPCA(A, options);
t_cpu = toc;
% calculate the adjusted variance explained by the extracted PCs
variance_rate = vars_adj(U,A);
fprintf('TPower: Propotion of explained variance: %f...\n', variance_rate);

% savw the results
TPower_Result = [];
TPower_Result.CPUTime = t_cpu;
TPower_Result.cardinality_vec = cardinality_vec;
TPower_Result.loadings = full(U);
TPower_Result.adjusted_variance_rate = variance_rate;

save ..\..\result\TPower_Result.mat 'TPower_Result';



