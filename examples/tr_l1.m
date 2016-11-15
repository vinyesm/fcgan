%
%        SYNTHETIC EXPERIMENT FOR REGULARIZED k-SPCA
%
% Regularized problem for simultaneaously sparse and low rank matrices
%
%       min_X |A-X|^2+lambda*|X|_1+mu*|X|_tr
%
%
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
p=100;
sigma=0;
lambda=1;
mu=1;

%% block diagonal covariance with 5 blocks of different sizes
block_sz=[30 20 15 10 10 5 5 5];
nblocs=length(block_sz);
a=cell(1,nblocs);
C=zeros(p);
deb=1;
for i=1:nblocs %nombre de blocs
    fin=deb+block_sz(i)-1;
    a{i}=sparse(p,1);
    a{i}(deb:fin,1)=randn(block_sz(i),1);
    a{i}=a{i}/norm(a{i});
    C=C+a{i}*a{i}';
    deb=deb+block_sz(i);
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
moy = zeros(1,p); % vector of means
Y = mvnrnd(moy, C+sigma^2*eye(p), n)';
S=cov(Y');

%%
%% fcgan
param.method='asqp';
param.lambda=lambda;
param.mu=mu;
param.max_nb_atoms=500; % Right now the code assumes that both are the same
param.max_nb_iter=500;
max_iter_fw=500;
param.epsStop=1e-5;
param.debug=false;
param.lmo=@lmo_chain_lgl;
param.ws=1;

%%
figure(1);
imagesc([abs(C(1:100,1:100)) ones(100,5) abs(S(1:100,1:100))]); colormap gray;
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
pbaspect([(2*100+5)/100 1 1]);

figure(2);
imagesc([(C(1:100,1:100)~=0) ones(100,5) (S(1:100,1:100)~=0)]); colormap gray;
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
pbaspect([(2*100+5)/100 1 1]);

%% algorithm
[x, as, hist, iter] = cgan_tr_l1(C,param);



% %% param
% 
% param.f=4;
% param.diag=0;
% param.PSD=true;
% param.max_nb_main_loop=20;
% param.powerIter=100;
% param.stPtPowerIter=1000;
% param.niterPS=10000;%5000
% param.epsStop=1e-8;
% param.PSdualityEpsilon=1e-3;
% param.k=0;
% param.PSmu=0; %strong convexity
% param.verbose=1;
% param.debug=0;
% param.sloppy=0;
% param.max_nb_atoms=param.max_nb_main_loop*param.niterPS;
% param.cardfun=inf*ones(1,p);
% param.cardfun(10)=1;
% %%
% 
% param.lambda=lambda;
% %% as quadprog
% param.opt='asqp';
% [Z_as, ActiveSet_as, hist_as, param_as, flaga_as, it_as] = cgan_spca(inputData,param);
% dg_as=hist_as.dg_sup;
% tt_as=hist_as.time_sup;
% fprintf('lambda=%f\n',lambda);
% fprintf('......tt=%f dg=%f\n',hist_as.time(end),hist_as.dg(end));
% 
% %% as quadprog
% param.sloppy=0;
% param.max_nb_main_loop=10;
% param.opt='proxbcd';
% param.nbMainLoop=20;
% [Z_bcd, ActiveSet_bcd, hist_bcd, param_bcd, flag_bcd,it_bcd] = cgan_spca(inputData,param);
% fprintf('......tt=%f dg=%f\n',hist_bcd.time(end),hist_bcd.dg(end));
% dg_bcd=hist_bcd.dg_sup;
% tt_bcd=hist_bcd.time_sup;
% 
% 
% %% Duality Gap Figure
% 
% xplot{1}=tt_bcd;
% xplot{2}=tt_as;
% yplot{1}=dg_bcd;
% yplot{2}=dg_as;
% colors = [0 0 0
%     1 0 0];
% legendStr={'bcd','cg'};
% 
% figure(1)
% for i=1:2
% loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i}); 
% hold on
% end
% title('duality gap lgl');
% legend('show','Location','southwest');
% grid on
% hold off
% 
% 
% xplot{1}=hist_bcd.time;
% xplot{2}=hist_as.time;
% yplot{1}=hist_bcd.dg;
% yplot{2}=hist_as.dg;
% colors = [0 0 0
%     1 0 0];
% legendStr={'bcd','cg'};
% 
% figure(2)
% for i=1:2
% loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i}); 
% hold on
% end
% title('duality gap lgl');
% legend('show','Location','southwest');
% grid on
% hold off
% 
