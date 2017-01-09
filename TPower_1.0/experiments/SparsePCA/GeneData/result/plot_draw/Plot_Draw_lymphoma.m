clear all;

figure;
clf; hold on

load('..\PathSPCA_Result_lymphoma.mat');
load('..\TPower_Result_lymphoma.mat');


xData = TPower_Result_lymphoma.cardinality_path;
yData{1} = TPower_Result_lymphoma.variance_rate;
yData{2} = PathSPCA_Result_lymphoma.variance_rate;
yData{3} = PathSPCA_Result_lymphoma.variance_rate_bound;

legendStr = {'TPower','PathSPCA','OptimalBound'};
prettyPlot(xData,yData,legendStr,'Sparse PCA','Cardinality','Propotion of Explained Variance',0, 'Lymphoma.eps');