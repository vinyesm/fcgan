function [x, f] = TPower_DkS(W, options, x0)

%%   Truncated power method for densest k-subgraph problem.
%
%     max x'*W*x    subject to x \in {0,1}^n, ||x||_0 <= k
%
% *** inputs ***
% - W:                      n x n graph weight matrix
% - options:                a structure stores user specified parameters which include:
%    -- verbose:            level of verbosity (0: no output, 1: final, 2: iter (default), 3: debug
%    -- cardinality:        target size of subgraph (default 10)
%    -- optTol:             optimality tolerance (default: 1e-6)
%    -- maxIter:            max number of iteration (default: 50)
%    -- initType:           initialization type (1: top-1 variance, 2: top-k variance (default))
% - x0:                     initialization vector
%
% *** outputs ***
% - x:                      n-dimensional subgraph indicator vector with k non-zeros
% - f:                      objective value at the output x
%
% Refer to:
%   Xiao-Tong Yuan, Tong Zhang, Truncated Power Method for Sparse Eigenvalue Problems, Technical Report, 2011
%
% Copyright (C) 2011/2012 - Xiao-Tong Yuan. 

%% Set Parameters
if nargin < 2
    options = [];
end
[verbose, cardinality, optTol, maxIter, initType] = ...
    myProcessOptions(...
    options, 'verbose', 2, 'cardinality', 10, 'optTol', 1e-8, 'maxIter', 50, 'initType', 1);

% Output Parameter Settings
if verbose >= 3
    fprintf('Running TPower Method for DkS...\n');
    fprintf('TPower Optimality Tolerance: %.2e\n', optTol);
    fprintf('TPower Maximum Number of Iterations: %d\n', maxIter);
end

%% Output Log
if verbose >= 2
    fprintf('TPower %10s %10s\n', ' Iteration ', ' Objective Val ');
end

%% Default initialization

if (nargin < 3)
    switch initType
        case 1
            [val,idx]=max(diag(W));
            x0 = zeros(size(W,1),1);
            x0(idx) = 1;
        case 2
            [val,idx]=sort(sum(W), 'descend');
            x0 = zeros(size(W,1),1);
            x0(idx(1:cardinality)) = 1;    end
end
    
x = sparse(x0);

% power step
s = W*x;
g = 2*s;
f = x'*s;

% truncate step
x_t = truncate_operator(g, cardinality);

f_old = f;

i = 1;

%% Main algorithmic loop
while i <= maxIter

    s_t = W*x_t;
    f_t = x_t'*s_t;
    f = f_t;
       
    % ensure the objective to be non-decreasing
    lambda = 1e-4;
    while(f < f_old - 1e-10)
        g_t = g + 2*lambda*x;
        x_t = truncate_operator(g_t, cardinality);
        s_t = W*x_t;
        f_t = x_t'*s_t;
        f = f_t;
        lambda = lambda * 10;
    end
    
    if (abs(f-f_old) < optTol )
        break;
    end
    
   
    % Output Log
    if verbose >= 2
        fprintf('TPower %10d %10f \n', i, f);
    end
    
    x = x_t;
    g = 2*s_t;
    x_t = truncate_operator(g, cardinality);
    
    f_old = f;
    i = i+1;
end
f = full(f);
end

%% Evaluate the truncate operator
function u = truncate_operator(v , k)

u = zeros(length(v), 1);
[val, idx] = sort(abs(v), 'descend');

u(idx(1:k)) = 1;

u = sparse(u);
end
