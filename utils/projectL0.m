function y = projectL0(x,k)
% fprintf('%d %d %d\n',size(x,1),size(x,2),k);
 [z ix] = sort(abs(x),'descend');
 y = zeros(size(x));
 y(ix(1:k)) = sparse(x(ix(1:k)));

end