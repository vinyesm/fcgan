%
%        SYNTHETIC EXPERIMENT FOR REGULARIZED k-SPCA
%
% Regularized problem for chain latent group Lasso 
%
%       min_x |Dx-y|+lambda*chain_lgl(x)
%
% Comparison with Block Coordinate Descent
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

%% Experimental Setting
% n: nb samples 
% p: dimension
% k: block size kxk
% lambda: regularization parameter
% sigma: noise
n=20;
p=150;
sigma=0;
lambda=.1;

%% covariance
nblocs=5;
overlap=5;
a=cell(1,nblocs);
C=zeros(p);
for i=1:nblocs %nombre de blocs
    a{i}=sparse(p,1);
    a{i}(((10-overlap)*(i-1)+1):((10-overlap)*(i-1)+10),1)=ones(10,1);
    a{i}=a{i}/norm(a{i});
    C=C+a{i}*a{i}';
end

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
inputData.Y=S;
inputData.X1=eye(p);
inputData.X2=eye(p);

%%
figure(3);
imagesc([abs(C(1:40,1:40)) ones(40,5) abs(S(1:40,1:40))]); colormap gray;
hTitle=title('   covariance                         empirical covariance');
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
pbaspect([(2*40+5)/40 1 1]);


%% param

param.f=4;
param.diag=0;
param.PSD=true;
param.max_nb_main_loop=20;
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
param.cardfun(10)=1;
%%

param.lambda=lambda;
%% as quadprog
param.opt='asqp';
[Z_as, ActiveSet_as, hist_as, param_as, flaga_as, it_as] = cgan_spca(inputData,param);
dg_as=hist_as.dg_sup;
tt_as=hist_as.time_sup;
fprintf('lambda=%f\n',lambda);
fprintf('......tt=%f dg=%f\n',hist_as.time(end),hist_as.dg(end));

%% as quadprog
param.sloppy=0;
param.max_nb_main_loop=10;
param.opt='proxbcd';
param.nbMainLoop=20;
[Z_bcd, ActiveSet_bcd, hist_bcd, param_bcd, flag_bcd,it_bcd] = cgan_spca(inputData,param);
fprintf('......tt=%f dg=%f\n',hist_bcd.time(end),hist_bcd.dg(end));
dg_bcd=hist_bcd.dg_sup;
tt_bcd=hist_bcd.time_sup;


%% Duality Gap Figure

xplot{1}=tt_bcd;
xplot{2}=tt_as;
yplot{1}=dg_bcd;
yplot{2}=dg_as;
colors = [0 0 0
    1 0 0];
legendStr={'bcd','cg'};

figure(1)
for i=1:2
loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i}); 
hold on
end
title('duality gap lgl');
legend('show','Location','southwest');
grid on
hold off


xplot{1}=hist_bcd.time;
xplot{2}=hist_as.time;
yplot{1}=hist_bcd.dg;
yplot{2}=hist_as.dg;
colors = [0 0 0
    1 0 0];
legendStr={'bcd','cg'};

figure(2)
for i=1:2
loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i}); 
hold on
end
title('duality gap lgl');
legend('show','Location','southwest');
grid on
hold off

