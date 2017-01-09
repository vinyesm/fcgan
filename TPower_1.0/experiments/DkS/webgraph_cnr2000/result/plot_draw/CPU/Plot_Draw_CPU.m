clear all;

figure;


load('..\..\TPower_Result.mat');
yData{1} = TPower_Result.cpu_time;
xData = TPower_Result.cardinality;

load('..\..\Greedy_Ravi_Result.mat');
yData{2} = Greedy_Ravi_Result.cpu_time;

load('..\..\Greedy_Feige_Result.mat');
yData{3} = Greedy_Feige_Result.cpu_time;

legendStr = {'TPower', 'Greedy-Ravi', 'Greedy-Feige'};
prettyPlot(xData,yData,legendStr,'Densest k-Subgraph','Cardinality','CPU Time (Sec.)',1, 'cpu.eps');


