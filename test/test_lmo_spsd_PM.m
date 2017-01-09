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
n=5;
A=randn(n);
A=(A+A')/2;
param.stPtPowerIter=1000;
param.powerIter=200;
param.cardfun=(1:(n)).^.2;

%%
param.PSD=true;
tic
[u_sym,k_sym,allVal_sym]=lmo_spsd(A,param);
tsym=toc;
display(['symmetric PM  done in ' num2str(tsym)]);
%param.stPtPowerIter=1;
%param.powerIter=1000;
param.PSD=false;
tic
[u,v,kuv, allVal_uv]=lmo_multi_spsd(A,param);
tmulti=toc;
display(['multi PM  done in ' num2str(tmulti)]);

%%
figure(1);clf;
plot(allVal_sym,'b', 'LineWidth',2);hold on
plot(allVal_uv,'r', 'LineWidth',2);
