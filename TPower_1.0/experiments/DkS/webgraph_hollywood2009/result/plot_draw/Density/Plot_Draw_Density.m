clear all;

figure;
clf; hold on

load('..\..\TPower_Result.mat');
yData{1} = TPower_Result.density;
xData = TPower_Result.cardinality;

load('..\..\Greedy_Ravi_Result.mat');
yData{2} = Greedy_Ravi_Result.density;

load('..\..\Greedy_Feige_Result.mat');
yData{3} = Greedy_Feige_Result.density;

legendStr = {'TPower', 'Greedy-Ravi', 'Greedy-Feige'};
prettyPlot(xData,yData,legendStr,'Densest k-Subgraph','Cardinality','Density',0, 'density.eps');


