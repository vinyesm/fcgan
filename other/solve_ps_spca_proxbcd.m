function [Z, D, ActiveSet, hist] = solve_ps_spca_proxbcd(Z,D,ActiveSet,param,inputData)


cont = 1;
count=0;
i = 1;

obj = zeros(1,param.niterPS);
pen =  zeros(1,param.niterPS);
loss =  zeros(1,param.niterPS);
dg =  zeros(1,param.niterPS);
time =  zeros(1,param.niterPS);

%precompute cardinality function values f(k)
cardVal=zeros(1, length(ActiveSet.k));
for ii=1:length(ActiveSet.tracenorm)
    cardVal(ii)=param.cardfun(ActiveSet.k{ii});
end

if param.diag==1
    %quantities used for optimization of diag
    S=inputData.X1*inputData.X2;
    S2=S.*S;
end
% YStart=inputData.Y;

while cont
    if isempty(ActiveSet.I)
        break;
    end
    
    
    count = count+1;
    %     fprintf( 'iteration %d, nb sets S=%d ----------------\n',i, length(ActiveSet.I));
    
    %% optimization wrt block ix
    deb=0;
    if count > deb % o choose if first optimizing D or block ix
        ix = randi(length(ActiveSet.I),1);
        currI = ActiveSet.I{ix};
        %         currJ = ActiveSet.J{ix};
        currJ=currI;
        currPassiveZ = Z;
        currPassiveZ(currI,currJ) = Z(currI,currJ) - ActiveSet.Z{ix};
        
        if param.f == 1       % prox : the subproblem indexed by ix is solved in one iteration, a singular value shrinkage
            if param.PSD
                [currZ, tracenorm, u, s] = ShrinkPSD(inputData.Y(currI,currI) - currPassiveZ(currI,currI) - diag(D(currI)) , param.lambda); v = u;
            else
                [currZ, tracenorm, u, v, s] = Shrink(inputData.Y(currI,currJ) - currPassiveZ(currI,currJ) - diag(D(currI)), param.lambda);
            end
            
        elseif param.f == 2   % multitask
            currY = inputData.Y - inputData.X*currPassiveZ;
            [currZ, u, v, s, grad, tracenorm] = fista(inputData.X(:,currI),currY(:,currJ), param, ActiveSet.Z{ix});
            
        elseif param.f == 3   % quadratic
            currY = inputData.Y - diag(inputData.X*currPassiveZ*inputData.X');
            [currZ, u, v, s, grad, tracenorm] = fista(inputData.X(:,currI) , currY , param, ActiveSet.Z{ix});
            
        elseif param.f == 4   % bilinear
            currY = inputData.Y - inputData.X1*currPassiveZ*inputData.X2 - inputData.X1*diag(D)*inputData.X2 ;
            X{1} = inputData.X1(:,currI);
            X{2} = inputData.X2(currJ,:);
            %             keyboard;
            [currZ, u, v, s, grad, tracenorm] = fista(X , currY , param, ActiveSet.Z{ix});
            %             keyboard;
            
        end
        
        ActiveSet.Z{ix} = currZ;
        ActiveSet.U{ix} = u;
        ActiveSet.V{ix} = v;
        ActiveSet.Sigma{ix} = s;
        Z = currPassiveZ;
        Z(currI,currJ) = Z(currI,currJ) + currZ;
        ActiveSet.tracenorm{ix} = tracenorm;
        ActiveSet.fronorm{ix} = norm(currZ,'fro')^2;
        
        if param.sloppy==0 || (param.sloppy==1 && mod(count,100)==1)
            [loss(i),pen(i),obj(i),dg(i),time(i)]=getvalProx(Z,D,ActiveSet,inputData,param,cardVal);            
            %% stopping criterion
            cont = (dg(i)>param.PSdualityEpsilon) && count< 2*param.niterPS;
            i=i+1;
        end
    end    
    count=count+1;
    %% optimization wr to diag
    if param.diag==1
        i=i+1;
        D=S2\diag(S-S*Z*S);
        D = max(D,0);
        
        if param.sloppy==0 || (param.sloppy~=0 && mod(count-1,100)==1)  
            [loss(i),pen(i),obj(i),dg(i),time(i)]=getvalProx(Z,D,ActiveSet,inputData,param,cardVal);
            %% stopping criterion
            cont = (dg(i)>param.PSdualityEpsilon) && count< 2*param.niterPS;
            i=i+1;
        end
    end
end

i=i-1;
%redimensioning arrays
hist.obj = obj(1:i);
hist.pen =  pen(1:i);
hist.loss =  loss(1:i);
hist.dg =  dg(1:i);
hist.time = time(1:i);
hist.obj_sup=obj(1);
hist.dg_sup=dg(1);
hist.time_sup=time(1);

% 
% i=i-1;
% obj = obj(1:i);
% pen =  pen(1:i);
% loss =  loss(1:i);
% dg =  dg(1:i);
% time = time(1:i);
% objSup = zeros(1,i);
% %dgSup = zeros(1,i);
% ttSup=time(1);
% objSup(1) = obj(1);
% dgSup(1) = dg(1);
% % nbit=count;


if param.debug==1 && i>0
    figure(17);clf
    subplot(1,2,1)
    idx=find(obj);
    plot(idx,obj(idx));
    subplot(1,2,2)
    semilogy(idx,dg(idx));
    fprintf('     solvePS i=%d  dg=%f/%f\n',i,dg(i),param.PSdualityEpsilon);
    keyboard;
end

if i>=param.niterPS
    fprintf('     solvePS i=%d  dg=%f/%f\n',i,dg(i),param.PSdualityEpsilon);
    
end

end



