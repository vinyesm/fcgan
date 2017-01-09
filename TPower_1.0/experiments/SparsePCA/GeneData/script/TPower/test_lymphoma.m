clear all;
load ..\..\data\LymphomaCov.mat

A = covlymph;

%%-------------------------------- TPower for Sparse PCA -------------------------------
% parameters setting
options = [];
options.optTol = 1e-4;
options.maxIter = 500;
options.verbose = 0;
options.initType = 1;

dim = size(A,1);
cardinality_path = 1:dim;
disp(['Computing Sparse Loadings with Cardinality Path: 1:', int2str(dim)]);

U = [];
vars = [];
tic; 
for ( k = cardinality_path )
    
    options.cardinality_vec = [k];
    [u f] = TPower_SPCA(A, options);
    U = [U, u];
    vars = [vars, f];
end
t_cpu = toc;
    
% save results
TPower_Result_lymphoma = [];
TPower_Result_lymphoma.CPUTime = t_cpu;
TPower_Result_lymphoma.variance_rate = vars/trace(A);
TPower_Result_lymphoma.cardinality_path = cardinality_path;
save ..\..\result\TPower_Result_lymphoma.mat 'TPower_Result_lymphoma';



