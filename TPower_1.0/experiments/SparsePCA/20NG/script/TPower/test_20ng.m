% Test code of TPower on 20NG
clear all; 
load ..\..\data\covariance20NGTop1K.mat

A = covng20_1K;

% Set up Objective Function

%% ------------------------------- TPower for Sparse PCA -------------------------------

% parameters setting
cardinality_vec = [20, 20, 10 10, 10];

options = [];
options.cardinality_vec = cardinality_vec;
options.OptTol = 1e-8;
options.MaxIter = 500;
options.verbose = 0;
options.initType = 1;

disp('TPower: Computing Sparse Loadings with Cardinality:');
disp(cardinality_vec);

% run TPower 
tic;
U = TPower_SPCA(A, options);
t_cpu = toc;
% calculate the adjusted variance explained by the extracted PCs
U = full(U);
variance_rate = vars_adj(U,A);
fprintf('TPower: Propotion of explained variance: %f...\n', variance_rate);

keywords = [];
for (i=1:size(U,2))
    u = U(:,i);
    keywords_cur = terms(find(u~=0));
    
    keywords = [keywords, {keywords_cur}];
end

% Save results
TPower_Result = [];
TPower_Result.CPUTime = t_cpu;
TPower_Result.cardinality_vec = cardinality_vec;
TPower_Result.variance_rate = variance_rate;
TPower_Result.keywords = keywords;


save ..\..\result\TPower_Result.mat 'TPower_Result';



