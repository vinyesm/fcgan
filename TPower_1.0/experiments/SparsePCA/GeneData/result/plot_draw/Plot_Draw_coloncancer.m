clear all;

figure;
clf; hold on

load('..\TPower_Result_coloncancer.mat');
load('..\PathSPCA_Result_coloncancer.mat');




xData = TPower_Result_coloncancer.cardinality_path;
yData{1} = TPower_Result_coloncancer.variance_rate;
yData{2} = PathSPCA_Result_coloncancer.variance_rate;
yData{3} = PathSPCA_Result_coloncancer.variance_rate_bound;

legendStr = {'TPower','PathSPCA','OptimalBound'};
prettyPlot(xData,yData,legendStr,'Sparse PCA','Cardinality','Propotion of Explained Variance',0, 'Colon_Cancer.eps');


