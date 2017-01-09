%%
%testing different cardinality functions

%% add paths
clc; clear all;
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');
addpath('../other');

%%
n=40+20+10+5;
A0=zeros(n);
W=randn(n);
W=(W+W')/2;
A0(1:40,1:40)=rand_sym_mat(40,40o);
A=A0+0*W;
param.stPtPowerIter=1000;
param.powerIter=200;

%%
alpha=[.2 .5 .7 1 2];

for i=1:length(alpha)
param.cardfun=(1:(n)).^alpha(i);
[uBest{i},kBest{i},allVal{i}]=lmo_spsd(A,param);
display(['alpha=' num2str(alpha(i)) '  done']);
end

%%
figure(1);clf;
for i=1:length(alpha)
plot(allVal{i},'Color',[0,0.2*(i-1),1],'LineWidth',2);hold on;
end