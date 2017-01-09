function [x] = Greedy_DkS_Feige(A, k)

%% Re-implementation of the greedy densest k-subgraph (DkS) algorithm as described in 
%
%    U. Feige, D. Peleg, G. Kortsarz, The Dense k-Subgraph Problem. Algorithmica 29(3): 410-421 (2001)
%
% By Xiao-Tong Yuan,  November 2011

% *** inputs ***
% - A:         n x n graph weight matrix
% - k          target size of subgraph 
%
% *** outputs ***
% - x:         n-dimensional  subgraph indicator vector with k non-zeros

%% The greedy procedure ( the procedure 2 in Feige's paper)
dim = size(A,1);

d = sum(A,2);

[val,idx] = sort(d, 'descend');

% Find the node set H with top k/2 degrees
s = floor(k/2);
H = idx(1:s);

% Find the node set C (from the remaining nodes) with top k/2 neiborhoods in H
d = sum(A(:, H),2);
d(H) = -inf;
[val,idx] = sort(d, 'descend');
C = idx(1:k-s);

active_set = [H;C];
x = zeros(dim,1);
x(active_set) = 1;
end
