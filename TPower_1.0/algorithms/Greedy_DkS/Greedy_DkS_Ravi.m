function [x] = Greedy_DkS_Ravi(A, k)

%% Re-implementation of the greedy densest k-subgraph (DkS) algorithm as described in 
%
%    S. S. Ravi, D. J. Rosenkrantz and G. K. Tayi, Heuristic and Special Case Algorithms for Dispersion Problems. Operations Research  42(2): 299-310  (1994)
%
% By Xiao-Tong Yuan,  November 2011

% *** inputs ***
% - A:         n x n graph weight matrix
% - k          target size of subgraph 
%
% *** outputs ***
% - x:         n-dimensional  subgraph indicator vector with k non-zeros


%% the algorithm procedure
dim = size(A,1);

if (dim < 1e6)
    [idx_x, idx_y, d] = find(A);
else
    [idx_x, idx_y, d] = find(A(1:1e4, 1:1e4));
end

[val,idx] = max(d);

active_set = [idx_x(idx), idx_y(idx)];
for (i=1:k-1)   

   d = sum(A(:, active_set),2);

   d(active_set) = -inf;
   [val, idx] = max(d);
   active_set = [active_set, idx];
   
end
x = zeros(dim,1);
x(active_set) = 1;
end
