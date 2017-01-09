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
A0(1:40,1:40)=rand_sym_mat(40,2);
A=A0+0*W;
param.stPtPowerIter=100;
param.powerIter=200;

%%
alpha=[0.01 0.05 0.1 0.5 0.8 1 2 5];

tic
for i=1:length(alpha)
param.cardfun=(1:(n)).^alpha(i);
[uBest{i},kBest{i},allVal{i}]=lmo_spsd_TPower(A,param);
display(['alpha=' num2str(alpha(i)) '  done']);
end
toc

% tic
% for i=1:length(alpha)
% param.cardfun=(1:(n)).^alpha(i);
% %[uBest{i},kBest{i},allVal{i}]=lmo_spsd_TPower(A,param);
% [uBest2{i},kBest2{i},allVal2{i}]=lmo_spsd(A,param);
% display(['alpha=' num2str(alpha(i)) '  done']);
% end
% toc

%%
figure(1);clf;
for i=1:length(alpha)
plot(allVal{i},'Color',[0,(i-1)/(length(alpha)-1+0.01),1],'LineWidth',2);hold on;
stem(kBest{i},allVal{i}(kBest{i}),'.','Color',[1,0,0],'LineWidth',2);hold on;
%plot(allVal2{i},'Color',[0.2*(i-1),0,1],'LineWidth',2);hold on;
end
