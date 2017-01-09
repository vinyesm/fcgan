function [loss,pen,obj,dualityGap,time]=get_val_lgm_asqp(Z,D,ActiveSet,inputData,param,cardVal)

YStart=inputData.Y;
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
