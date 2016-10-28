%
%        SYNTHETIC EXPERIMENT FOR REGULARIZED k-CHAIN LASSO
%
% Regularized problem for chain latent group Lasso 
%
%       min_x |Dx-y|+lambda*chain_lgl(x)
%
% Comparison with Frank-Wolfe
%
% 
%%%%%%%%%%%%
% Marina Vinyes and Guillaume Obozinski, 2016
% %%%%%%%%%%%%

%% add paths
clc; clear all;
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');

%% Experimental Setting
% n: nb samples 
% p: dimension
% k: chain length
% lambda: regularization parameter
% sigma: noise
n=300;
p=1000;
k=8;
lambda=5;
sigma=0.1;

%% synthetic data
% Random design matrix D, paramter vector x, output y
m=10;
x=[2*ones(m,1);zeros(p-m,1)];

U=randn(n,n);
[U,R]=qr(U);
V=randn(p,p);
[V,R]=qr(V);
S=zeros(n,p);
s=0.95.^(0:1:(min(n,p)-1));
S(1:min(n,p),1:min(n,p))=diag(s);
D=sqrt(n)*U*S*V';
 
y=D*x+sigma*randn(n,1);

%% setting
% k=7;
%param.x0=x;
%param.debug=true;
param.lmo=@lmo_chain_lgl;
param.flag_no_design=false;

%% fcgan

param.lambda=lambda;
param.K=k;
param.max_nb_atoms=500; % Right now the code assumes that both are the same
param.max_nb_iter=500;
max_iter_fw=500;
param.epsStop=1e-5;
param.debug=false;
param.lmo=@lmo_chain_lgl;

%% experiments

% Column Generation Algoritm
param.method='asqp';
tic
[x_sol, ActiveSet, hist, iter] = cgan_lgl(D,y,param);
cg_x=x_sol;
cg_dg=hist.dg;
cg_tt=hist.tt;
fprintf(['colgen solution is ' num2str(sum(x_sol~=0)) '-sparse\n'])


% Frank-Wolfe with decrasing stepsize
param.method='FW';
tic
param.max_nb_iter=max_iter_fw;
param.max_nb_atoms=max_iter_fw;
[x_sol, ActiveSet, hist, iter] = cgan_lgl(D,y,param);
fw_x=x_sol;
fw_dg=hist.dg;
fw_tt=hist.tt;
fprintf(['FW solution is ' num2str(sum(x_sol~=0)) '-sparse\n'])

% Line Search Frank-Wolfe
param.method='FWls';
tic
param.max_nb_iter=max_iter_fw;
param.max_nb_atoms=max_iter_fw;
[x_sol, ActiveSet, hist, iter] = cgan_lgl(D,y,param);
fwls_x=x_sol;
fwls_dg=hist.dg;
fwls_tt=hist.tt;
fprintf(['FWls solution is ' num2str(sum(x_sol~=0)) '-sparse\n'])

% Pair Wise Frank-Wolfe
param.method='FWpw';
tic
param.max_nb_iter=max_iter_fw;
param.max_nb_atoms=max_iter_fw;
[x_sol, ActiveSet, hist, iter]  = cgan_lgl(D,y,param);
fwpw_x=x_sol;
fwpw_dg=hist.dg;
fwpw_tt=hist.tt;
fprintf(['FWpw solution is ' num2str(sum(x_sol~=0)) '-sparse\n'])





%% Duality Gap Figure

xplot{1}=fw_tt;
xplot{2}=fwls_tt;
xplot{3}=fwpw_tt;
xplot{4}=cg_tt;
yplot{1}=fw_dg;
yplot{2}=fwls_dg;
yplot{3}=fwpw_dg;
yplot{4}=cg_dg;
colors = [1 0 1
    .5 1 .5    
    0 1 1
    1 0 0];
legendStr={'FW', 'FW-ls', 'FW-pw','colgen'};

figure(1)
for i=1:4
loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i}); 
hold on
end
title('duality gap lgl');
legend('show','Location','southwest');
grid on
hold off
