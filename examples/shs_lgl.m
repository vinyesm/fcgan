%EXPERIMENT WHS, SINTHETIC DATA
%
% weak hierachical sparsity wij~=0 => wi~=0 ou wj~=0
% Lasso formulation from Hierarchical sparse modeling: A choice of two regularizers 
% Yan, X. and Bien, J. (2015) - arXiv preprint arXiv:1512.01631
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
% p: dimension, effective dimension is p*(p-1)/2+p
% s: sparsity parameter
% nb_interactions
% lambda: regularization parameter
% sigma: noise
n=1000;
p=50; 
s=0.1;
nb_interactions=ceil(s*p*(p-1)/2); % < p*p(-1)/2
sigma=0.1;
lambda=1*n;

%groups defining the norm whs
%groups
groups=cell(p+(p*(p-1))/2,1);
iter=0;
for i=1:p
    iter=iter+1;
    groups{iter}=i;
end
for i=1:p
    for j=(i+1):p
        iter=iter+1;
        groups{iter}=[i,j];
    end
end
%group_mat, use in lmo (selection of the new atom)
group_mat=sparse(p+p*(p-1)/2,p+p*(p-1));
gg=cell2mat(groups((p+1):end));
group_mat(1:p,1:p)=speye(p);
index=p+1;
for i=1:p
    J=index:(index+p-2);
    group_mat(i,J)=1;
    ii=sum(gg==i,2);
    JJ=find(ii);
    group_mat((p+1):end,J)=sparse(JJ,1:(p-1),1,p*(p-1)/2,p-1);
    index=index+p-1;
end

%% synthetic data
% random Design matrix X_all, paramter vector w0, output Y_all

pg=length(groups);
group2=cell2mat(groups(p+1:end));
idx=[ones(1,nb_interactions),zeros(1,pg-p-nb_interactions)];
idx=idx(randperm(pg-p));

w0=sparse(pg,1);
w0(idx>0)=1;
w0(unique(group2(idx>0,2)))=0.5;

X=randn(n,p);
X_all=X;
for i=1:(p-1)
    X_all=[X_all bsxfun(@times,X(:,i),X(:,(i+1):end))];
end

Y_all=X_all*w0+sigma*randn(n,1);

%% fcgan

param.lambda=lambda;
param.K=2;
param.max_nb_atoms=5000; 
param.max_nb_iter=5000;
param.epsStop=1e-5;
param.debug=false;
param.lmo=@lmo_whs;
param.flag_no_design=0;
param.group_mat=group_mat;
param.p0=p;
param.method='asqp';

%% Experiments
fprintf(' whs n=%d p=%d\n',n,p);

%with Warm Start in asqp
param.p0=p;
param.ws=1;

[x, ActiveSet, hist] = cgan_lgl(X_all,Y_all,param);
x_as=x;
dg_as=hist.dg;
tt_as=hist.tt;
fprintf('.........asqp  tt=%d  dg=%f\n',tt_as(end),dg_as(end));

%without Warm Start in asqp
param.p0=p;
param.ws=0;
[x, ActiveSet, hist] = cgan_lgl(X_all,Y_all,param);
x_as0=x;
dg_as0=hist.dg;
tt_as0=hist.tt;
fprintf('.........asqp no ws  tt=%d  dg=%f\n',tt_as0(end),dg_as0(end));

%% Duality Gap Figure

xplot{1}=tt_as0;
xplot{2}=tt_as;
yplot{1}=dg_as0;
yplot{2}=dg_as;
colors = [0 0 0
    1 0 0];
legendStr={'colgen without ws','colgen with ws'};

figure(1)
for i=1:2
    loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i});
    hold on
end
title('duality gap whs regularized');
legend('show','Location','southwest');
grid on
hold off
