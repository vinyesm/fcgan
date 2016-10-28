function [loss,pen,obj,dualityGap,time]=get_val_spca_proxbcd(Z,D,ActiveSet,inputData,param,cardVal)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
YStart=inputData.Y;

p=size(Z,1);
if param.f == 1
    currloss =  .5*norm(Z - inputData.Y, 'fro')^2;
elseif param.f == 2
    currloss =  .5*norm(inputData.X*Z - inputData.Y, 'fro')^2;
elseif param.f == 3
    currloss = .5*norm(inputData.Y - diag(inputData.X*Z*inputData.X'))^2;
elseif param.f == 4
    currloss = .5*norm(inputData.Y - inputData.X1*(Z+diag(D))*inputData.X2, 'fro')^2;
end
loss = currloss;
pen = sum(cell2mat(ActiveSet.tracenorm).*cardVal);
obj = currloss + param.PSmu/2*sum(cell2mat(ActiveSet.fronorm)) + param.lambda*pen;
inputData.Y= YStart - inputData.X1*diag(D)*inputData.X2;
dualityGap = get_dg_spca_asqp(Z,param,ActiveSet,inputData,obj,loss,pen);
inputData.Y=YStart;
time=toc;
end

