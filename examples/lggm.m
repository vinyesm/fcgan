%
%        SYNTHETIC EXPERIMENT FOR LEARNING GGM WITH LATENT VARIABLES
%
% Regularized problem for chain latent group Lasso 
%
%       min_x |S^.5*Z*S^.5-S|+lambda*omega(Z)
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
addpath('../other');
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');

%% Experimental Setting
n=500;
p=15;
sigma=0;
lambda=0.5;

%% Covariance and design matrix
rho1 = 0.7;
rho2 = 0.9;
C=eye(15);
C(1:5,1:5) = rho1*ones(5)+(1-rho1)*eye(5);
C(6:10,6:10) = rho2*ones(5)+(1-rho2)*eye(5);
des=inv(C);

%% output
out_as.cpu=[];
out_as.na=[];
out_as.qp=[];
out_as.hessian=[];
out_as.time=[];
out_bcd.cpu=[];
out_cndg.cpu=[];

%% data
mu = zeros(1,p); % vector of means
Y = mvnrnd(mu, C+sigma^2*eye(p), n)';
S=cov(Y');

%%
figure(1);
imagesc([abs(C) 1*ones(p,5) abs(inv(S))]); colormap gray;
hTitle=title('   inverse covariance                        inverse empirical covariance');
% set(gca,'XTick',0:20:135,'YTick',0:10:40);
h1=xlabel('');
h2=ylabel('');
set(h1,'Visible','off');
set(gca,'FontName','AvantGarde','FontWeight','normal','FontSize',12);
set(hTitle,'FontName','AvantGarde','FontSize',14,'FontWeight','bold');
set(h2,'FontName','AvantGarde','FontSize',14,'FontWeight','normal');
set(gca,'xtick',[])
caxis([0, 0.3]);
colorbar;
pbaspect([(2*p+5)/p 1 1]);

keyboard;
%% param

param.f=4;
param.diag=0;
param.PSD=true;
param.max_nb_main_loop=100;
param.powerIter=100;
param.stPtPowerIter=1000;
param.niterPS=10000;%5000
param.epsStop=1e-8;
param.PSdualityEpsilon=1e-3;
param.k=0;
param.PSmu=0; %strong convexity
param.verbose=1;
param.debug=0;
param.sloppy=0;
param.max_nb_atoms=param.max_nb_main_loop*param.niterPS;
param.cardfun=inf*ones(1,p);
param.cardfun(5)=1;
%%
%S=C;
inputData.X1=S^.5;
inputData.X2=inputData.X1;
inputData.Y=eye(p);
param.lambda=lambda;
%% as quadprog
param.opt='asqp';
[Z_as,D_as, ActiveSet_as, hist_as, param_as, flaga_as, it_as] = cgan_lgm(inputData,param);
dg_as=hist_as.dg_sup;
tt_as=hist_as.time_sup;
fprintf('lambda=%f\n',lambda);
fprintf('......tt=%f dg=%f\n',hist_as.time(end),hist_as.dg(end));


%% Duality Gap Figure

xplot{1}=tt_as;
yplot{1}=dg_as;
colors = [1 0 0];
legendStr={'cg'};

figure(2);clf
for i=1
loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i}); 
hold on
end
title('duality gap lgl');
legend('show','Location','southwest');
grid on
hold off



xplot{1}=hist_as.time;
yplot{1}=hist_as.dg;
colors = [    1 0 0];
legendStr={'cg'};

figure(3);clf
for i=1
loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i}); 
hold on
end
title('duality gap lgl');
legend('show','Location','southwest');
grid on
hold off

figure(4);clf
loglog(tt_as,dg_as,'-','LineWidth',2,'Color',[0 0 0],'DisplayName','dg sup');hold on;
loglog(hist_as.time,hist_as.dg,'-','LineWidth',2,'Color',[1 0 0],'DisplayName','dg sup');hold on;
legend('show','Location','southwest');
grid on
hold off

%%
%%
figure(5);
imagesc([abs(C) 1*ones(p,5) abs(inv(S)) 1*ones(p,5) abs(inv(Z_as+diag(D_as)))]); colormap gray;
hTitle=title('   inverse covariance                        inverse empirical covariance      estimated inverse');
% set(gca,'XTick',0:20:135,'YTick',0:10:40);
h1=xlabel('');
h2=ylabel('');
set(h1,'Visible','off');
set(gca,'FontName','AvantGarde','FontWeight','normal','FontSize',12);
set(hTitle,'FontName','AvantGarde','FontSize',14,'FontWeight','bold');
set(h2,'FontName','AvantGarde','FontSize',14,'FontWeight','normal');
set(gca,'xtick',[])
caxis([0, 0.3]);
colorbar;
pbaspect([(3*p+10)/p 1 1]);

%%
figure(6);clf
nplots=length(ActiveSet_as.matrix_atoms)+2;
ncol=ceil(sqrt(nplots));
nrow=ceil(nplots/ncol);
iplot=1;
M=zeros(p);
while iplot<nplots-1
subplot(nrow,ncol,iplot);
imagesc(abs(ActiveSet_as.alpha(iplot)*ActiveSet_as.matrix_atoms{iplot}));
M=M+ActiveSet_as.alpha(iplot)*ActiveSet_as.matrix_atoms{iplot};
iplot=iplot+1;
end
subplot(nrow,ncol,iplot);
imagesc(abs(M));
subplot(nrow,ncol,iplot);
imagesc(abs(M+diag(D_as)));