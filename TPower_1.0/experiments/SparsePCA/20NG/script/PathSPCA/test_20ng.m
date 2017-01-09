% Test code of PathSPCA on 20NG
clear all;

load ..\..\data\covariance20NGTop1K.mat
S=covng20_1K;

cardinality_vec = [20, 20, 10 10, 10];
disp('PathSPCA: Computing Sparse Loadings with Cardinality:');
disp(cardinality_vec);
% run GPower 
tic;
[V]= PartialPathCov_Multiple(S, cardinality_vec);%FullPathCov(S);
t_pathspca = toc;
% [bnds,rhov,dualvals]=UpperBounds(F,S,res(:,1:kp));


% % Plot with optimal points. Uncomment to merge multiple plots
% epsr=1e-4;
% error=max(bnds(1:kp)-vars(1:kp),zeros(size(vars(1:kp))));
% optpi=find((error./vars(1:kp))<=epsr);optpv=vars(find((error./vars(1:kp))<=epsr));
% %errorl=max(bndsl-varsl(1:kp),zeros(size(vars(1:kp))));
% %optpil=find((errorl./varsl(1:kp))<=epsr);optpvl=varsl(find((errorl./varsl(1:kp))<=epsr));
% plot(1:kp,vars(1:kp),'-b',1:kp,bnds,':k','LineWidth',1.5);hold on;
% %plot(1:kp,varsl(1:kp),'-b',1:kp,bndsl,':k','LineWidth',1.5);
% %plot(optpil,optpvl,'b.','MarkerSize',20);
% plot(optpi,optpv,'b.','MarkerSize',20);
% xlabel('card');ylabel('var');
% %semilogy(error./vars,'bo');xlabel('Card');ylabel('Relative Error');axis([1 kp 1e-6 1e2]);

variance_rate = vars_adj(V, S);

fprintf('PathSPCA: Propotion of explained variance: %f...\n', variance_rate);
% bnds_rate = bnds / trace(covcolon);

keywords = [];
for (i=1:size(V,2))
    v = V(:,i);

    keywords_cur = terms(find(v~=0));
    
    keywords = [keywords, {keywords_cur}];
end

% Save results
PathSPCA_Result = [];
PathSPCA_Result.CPUTime = t_pathspca;
PathSPCA_Result.cardinality_vec = cardinality_vec;
PathSPCA_Result.variance_rate = variance_rate;
PathSPCA_Result.keywords = keywords;

save ..\..\result\PathSPCA_Result.mat 'PathSPCA_Result';