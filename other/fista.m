function [S u v s grad tracenorm] = fista(X,Y, param,Start)
if param.f < 4
    L = 2*norm(X'*X);
else 
    L = 2*norm(X{1})^2*norm(X{2})^2 ;
end
S = Start;
Sx = Start;
tk = 1;
tkp1 = tk;

 if param.f == 3
        L = L^2;
 end

for iter = 1:param.innerLoopIter
    
    if param.f == 2
        grad = X'*(X*S-Y);
    elseif param.f ==3
        grad = X'*diag(diag(X*S*X')-Y)*X;
    elseif param.f == 4
%         grad = X{1}'*(X{1}*S*X{2}-Y)*X{2}';
          grad = X{1}'*(X{1}*S*X{2}-Y)*X{2}' + param.PSmu*S;
    end
        
        SxOld = Sx;     
        Sx = S - 1/L*grad;

        if param.PSD
         card = size(X{1},2);
%          fprintf('card=%d\n',card);
%          keyboard;
         [Sx tracenorm u s] = ShrinkPSD(Sx, param.lambda*param.cardfun(card)  / L ); v = u;
%          keyboard;
%          [Sx tracenorm u s] = ShrinkPSD(Sx, param.lambda  / L ); v = u;
        else
          [Sx tracenorm, u v s] = Shrink(Sx,param.lambda / L);
        end 
        tk = tkp1;
        tkp1 = .5*(1+sqrt(1+4*tk^2));
        
        S = Sx + (tk-1)/tkp1*(Sx-SxOld);
        
        % only a prox step
        S = Sx; 
        break;
end

end