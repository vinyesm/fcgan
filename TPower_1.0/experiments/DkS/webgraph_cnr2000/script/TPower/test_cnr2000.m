clear all;
load ..\..\data\data.mat


%% -------------------------------- TPower for finding DkS with different size k -------------------------------
cadinality = [100 500 1000:1000:10000];
W = ( W + W' )/2;

density_all = [];
cpu_time_all = [];
% Output Log
fprintf('TPower %10s %10s      %10s \n', ' Cardinality (k) ', ' Density ', ' CPU Time ');

% Finding DkS with different size k
for (k = cadinality)    
       
    options = [];
    options.OptTol = 1e-6;
    options.cardinality = k;
    options.MaxIter = 500;
    options.verbose = 0;
    options.initType = 2;  % initialize with the top k degreed nodes
    
    dim = size(W,1);
    
    tic;
    [u, d] = TPower_DkS(W, options);
    t_cpu = toc;
    
    Density = d / k;    
    fprintf('TPower %10d        %10f     %10f  \n', k, Density, t_cpu);
    density_all = [density_all, Density];
    cpu_time_all = [cpu_time_all, t_cpu];
end

% Save results
TPower_Result = [];
TPower_Result.cpu_time = cpu_time_all;
TPower_Result.cardinality = cadinality;
TPower_Result.density = density_all;


save ..\..\result\TPower_Result.mat 'TPower_Result';


