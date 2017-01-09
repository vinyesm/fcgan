clear all; 
load ..\..\data\data.mat

%% -------------------------------- Greedy-Feige for finding DkS with different size k -------------------------------
cardinality = [10:10:300];

W = double(W);

density_all = [];
cpu_time_all = [];

% Output Log
fprintf('Greedy-Feigei %10s %10s      %10s \n', ' Cardinality (k) ', ' Density ', ' CPU Time ');

% Finding DkS with different size k
for (k = cardinality)
    tic;
    u = Greedy_DkS_Feige(W, k);
    t_cpu = toc;
    Density = u'*W*u / k;
    fprintf('Greedy-Feige %10d        %10f     %10f  \n', k, Density, t_cpu);
    density_all = [density_all, Density];
    cpu_time_all = [cpu_time_all, t_cpu];
end

% Save results
Greedy_Feige_Result = [];
Greedy_Feige_Result.cpu_time = cpu_time_all;
Greedy_Feige_Result.cadinality = cardinality;
Greedy_Feige_Result.density = density_all;

save ..\..\result\Greedy_Feige_Result.mat 'Greedy_Feige_Result';
   



