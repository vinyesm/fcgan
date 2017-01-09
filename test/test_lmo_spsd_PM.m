%%
%testing different cardinality functions

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

%%
n=10;
A0=randn(n);
A0=(A0+A0')/2;
sa=eigs(A0,1,'sa');
if sa<0
    A=A0-2*sa*eye(n);
else
    A=A0;
end
%keyboard;

param.stPtPowerIter=100;
param.powerIter=200;
param.cardfun=(1:n).^.1;

%%
param.PSD=true;
tic
[u_sym,k_sym,allVal_sym]=lmo_spsd(A,param);
tsym=toc;
display(['symmetric PM  done in ' num2str(tsym) ' val ' num2str(allVal_sym(k_sym)) ]);

% param.stPtPowerIter=1;
%param.powerIter=1000;
param.PSD=false;
tic
[u,v,kuv, allVal_uv]=lmo_multi_spsd(A,param);
tmulti=toc;
display(['multi PM  done in ' num2str(tmulti) ' val ' num2str(allVal_uv(kuv)) ]);


tic
[u_tp,ktp,allVal_tp]=lmo_spsd_TPower(A,param);
ttp=toc;
display(['TPower done in ' num2str(ttp) ' val ' num2str(allVal_tp(ktp)) ]);

%%
% %Tpower
% options.verbose=0;
% options.optTol=1e-8;
% options.maxIter=1000;
% options.cardinality_vec=n:-1:1;
% options.initType=1;
% %options.cardfun=param.cardfun;
% tic
% [X, F] = TPower_SPCA(A, options);
% ttpowerspca=toc;
% display(['TPower_SPCA  done in ' num2str(ttpowerspca)]);
%%
figure(1);clf;
plot(allVal_sym,'b', 'LineWidth',2);hold on
plot(allVal_uv,'r', 'LineWidth',2);
plot(allVal_tp,'g', 'LineWidth',2);

%[allVal_uv allVal_sym F']
