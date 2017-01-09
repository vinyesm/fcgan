clear all;
load ..\..\data\data.mat


%% -------------------------------- Greedy-Ravi for finding DkS with different size k -------------------------------

cardinality = [10:10:300];

W = double(W);

density_all = [];
cpu_time_all = [];

% Output Log
fprintf('Greedy-Ravi %10s %10s      %10s \n', ' Cardinality (k) ', ' Density ', ' CPU Time ');

for (k = cardinality)
    tic;
    u = Greedy_DkS_Ravi(W, k);
    t_cpu = toc;
    Density = u'*W*u / k;
    fprintf('Greedy-Ravi %10d        %10f     %10f  \n', k, Density, t_cpu);
    density_all = [density_all, Density];
    cpu_time_all = [cpu_time_all, t_cpu];
end

Greedy_Ravi_Result = [];
Greedy_Ravi_Result.cpu_time = cpu_time_all;
Greedy_Ravi_Result.cadinality = cardinality;
Greedy_Ravi_Result.density = density_all;

save ..\..\result\Greedy_Ravi_Result.mat 'Greedy_Ravi_Result';






