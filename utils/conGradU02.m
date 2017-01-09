function x = conGradU02(A,nit,k)
% also named GPower or TPower
x = randn(size(A,1),1);
x = projectL0(x,k);
x = sparse(x);
debug=0;
% vold=0;
% vnew=0;
% x0 = zeros(size(A,1),1);
% if k==1
%     [val,idx]=max(diag(A));
%     x0(idx) = 1;
% else
%     [val,idx]=sort(diag(A), 'descend');
%     x0(idx(1:k)) = 1;
%     x0 = x0 / norm(x0);
% end
% x=x0;

if A ==0
    x = projectL0(x,k);
else
    
    for i=1:nit
        %             if test==1
        vold=x'*A*x;
        %             end
        x = A*x;
        x = projectL0(x,k);
        x = x/norm(x);
        %             if test==1
        vnew=x'*A*x;
        err=abs(vnew-vold)/vold;
        if (err)<1e-8
            if debug
                %                         keyboard;
                fprintf('        congrad converged in %d iteration\n', i);
            end
            break;
            %                 end
        end
    end
    if i>=nit
        %             keyboard;
        if debug
            fprintf('        congrad NOT converged at 0.000001. Only %f\n', err);
        end
    end
end
end
