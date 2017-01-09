clear all;
load ..\..\data\data.mat

%% -------------------------------- TPower for finding single DkS -------------------------------
W = ( W + W' )/2;
dim = size(W,1);

options = [];
options.OptTol = 1e-8;
options.MaxIter = 100;
options.verbose = 0;
options.initType = 2;  % initialize with the top k degreed nodes
    
cardinality_sequence = 30;

idx = 1:dim;
W_cur = W;
dense_subgraph_idx_mul = [];
center_city_idx_mul = [];


% Output Log
fprintf('TPower %10s %10s %10s     %10s \n', ' Subgraph ID ', ' Cardinality ', ' Density ', ' CPU Time ');

for (i = 1:length(cardinality_sequence))
    
    k = cardinality_sequence(i);
    options.cardinality  = k;
      
    
    tic;
    [u, d] = TPower_DkS(W_cur, options); %Power_Sparse_Proximal(W_cur, options,u0);
    t_cpu = toc;
    
    Density = d / k;     

    fprintf('TPower %10d   %10d     %10f     %10f  \n', i, k, Density, t_cpu);
    
    idx_1 = idx(find(u~=0));
    idx = idx(find(u==0));
    W_cur = W(idx, idx);
    
    city_name(idx_1)';
    
    deg = sum(W(idx_1, idx_1),2);
    [val , id] = max(deg);
    center_city_idx_mul= [center_city_idx_mul, idx_1(id)];
    
    dense_subgraph_idx_mul = [dense_subgraph_idx_mul, {idx_1}];
    
end

% Save results
TPower_Single_Result = [];
TPower_Single_Result.cardinality = cardinality_sequence;
TPower_Single_Result.dense_subgraph_idx = dense_subgraph_idx_mul;
TPower_Single_Result.center_city_idx = center_city_idx_mul;

save ..\..\result\TPower_Single_Result.mat 'TPower_Single_Result';
