function [loss,pen,obj,dualityGap,time]=get_val_spca_asqp(ActiveSet,inputData,param,cardVal)

YStart=inputData.Y;

p=size(YStart,1);

alph=ActiveSet.alpha;
smallValues=find(alph<0); % for numerical issues
alph(smallValues)=zeros(length(smallValues),1);    
ActiveSet.alpha=sparse(alph);

%% Update ActiveSet and Z and D    
Z=zeros(p);
nz=find(ActiveSet.alpha>1e-15);
for j=nz'
    u=ActiveSet.atoms(:,j);
    Z=Z+ActiveSet.alpha(j)*(u*u');
end

if param.diag==1
    D=ActiveSet.alpha(1:p);
else
    D=zeros(1,p);
end
Z=Z-diag(D);
inputData.Y= YStart - inputData.X1*diag(D)*inputData.X2;
currloss = .5*norm(inputData.Y - inputData.X1*Z*inputData.X2, 'fro')^2;
loss= currloss;
pen= dot(ActiveSet.alpha,cardVal);
obj = currloss + param.PSmu/2*sum(cell2mat(ActiveSet.fronorm)) + param.lambda*pen ;
dualityGap = get_dg_spca_asqp(Z,param,ActiveSet,inputData,obj,loss,pen);
dualityGap= real(dualityGap);
time=toc;
inputData.Y= YStart;

end

