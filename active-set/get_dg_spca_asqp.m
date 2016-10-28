function [dualityGap] = get_dg_spca_asqp(Z,param,ActiveSet,inputData,obj,loss,pen)

H = gradient(Z,inputData,param); % gradient
temp = -1;
shrink=-1;

if param.PSD
    
    if param.PSmu>0
        
        psistar=0;
        for i = 1:length(ActiveSet.I)
            currS=-H(ActiveSet.I{i}, ActiveSet.I{i});
            currS= currS-param.lambda*param.cardfun(ActiveSet.k{i})*eye(length(ActiveSet.I{i}));
            [Sh,tracenorm,u,s] = ShrinkPSD(currS,0);
            Sh=1/param.PSmu*Sh;
            psistar = psistar + trace(Sh'*currS) - param.PSmu/2*norm(Sh,'fro')^2;
        end
        
    else
        
        for i = 1:length(ActiveSet.I)
            currS=-H(ActiveSet.I{i}, ActiveSet.I{i});
            currS=0.5*(currS+currS');
            %         currTemp = PowerIteration(currS,param);
            if rcond(currS)<1e-14
%                 fprintf('bad conditioned\n');
            end
            if sum(sum(isnan(currS)))
                error('nan\n');
            end
            currTemp = eigs(currS,1,'la');%max(eig(currS));
            currTemp = currTemp/(param.cardfun(ActiveSet.k{i}));
            %         currTemp =  norm(H(ActiveSet.I{i}, ActiveSet.I{i})/param.cardfun(ActiveSet.k{i});
            
            if currTemp>temp
                temp = currTemp;
            end
            
        end
    end
    
else
    
    for i = 1:length(ActiveSet.I)
        currTemp = norm(H(ActiveSet.I{i}, ActiveSet.J{i}));
        if currTemp>temp
            temp = currTemp;
        end
    end
end

%test if temp<0 constraint already satisfied. We want to select 1 in min(1, param.lambda / temp)
if temp<0
    shrink=1;
else
    shrink=min(1, param.lambda / temp);
end

%duality gap
if param.PSmu>0
    kappa = inputData.X1*Z*inputData.X2-inputData.Y;
    dualityGap = obj + sum(sum( inputData.Y .* kappa )) + 1 / 2 * sum(sum( kappa .* kappa )) + psistar;
    gapLoss = loss + sum(sum( inputData.Y .* kappa )) + 1 / 2 * sum(sum( kappa .* kappa )) -sum(sum((inputData.X1*Z*inputData.X2).*kappa));
    gapPen = param.lambda*pen + psistar + sum(sum((inputData.X1*Z*inputData.X2).*kappa));
else
    kappa = shrink * (inputData.X1*Z*inputData.X2-inputData.Y);
    dualityGap = obj + sum(sum( inputData.Y .* kappa )) + 1 / 2 * sum(sum( kappa .* kappa ));
    gapLoss = loss + sum(sum( inputData.Y .* kappa )) + 1 / 2 * sum(sum( kappa .* kappa )) -sum(sum((inputData.X1*Z*inputData.X2).*kappa));
    gapPen = param.lambda*pen + sum(sum((inputData.X1*Z*inputData.X2).*kappa));
    
end

if dualityGap<0 && abs(dualityGap)>1e-10
        fprintf('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
%     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    dualityGap=abs(dualityGap);
end

end

