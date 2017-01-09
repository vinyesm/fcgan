function [X, F] = TPower_SPCA(A, options)

%%   Truncated power method for sparse principal component analysis problem.
%
%     max x'*A*x    subject to ||x||=1, ||x||_0 <= k
%
% *** inputs ***
% - A:                      p x p covariance matrix
% - options:                a structure stores user specified parameters which include:
%    -- verbose:            level of verbosity (0: no output, 1: final, 2: iter (default), 3: debug
%    -- cardinality_vec:    an m-dimensional vector stores the cadinality for each sparse loading (default [10])
%    -- optTol:             optimality tolerance (default: 1e-6)
%    -- maxIter:            max number of iteration (default: 50)
%    -- initType:           initialization type (1: top-1 variance, 2: top-k variance (default))
% 
% *** outputs ***
% - x:                      p-dimensional sparse vector with k non-zeros
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
[verbose, cardinality_vec, optTol, maxIter, initType] = ...
    myProcessOptions(...
    options, 'verbose', 2, 'cardinality_vec', 10, 'optTol', 1e-6, 'maxIter', 50, 'initType', 1);

% Output Parameter Settings
if verbose >= 3
    fprintf('Running TPower Method for Sparse PCA...\n');
    fprintf('TPower Optimality Tolerance: %.2e\n', optTol);
    fprintf('TPower Maximum Number of Iterations: %d\n', maxIter);
end

%% Output Log
if verbose >= 2
    fprintf('TPower %10s %10s\n', ' Iteration ', ' Objective Val ');
end

%% Initialization
X = [];
F = [];
dim = size(A,1);
m = length(cardinality_vec);

%% Main loop to extract multiple sparse loadings
for ( l = 1:m )
    
    cardinality = cardinality_vec(l);
    options_cur = [];
    options_cur.cardinality = cardinality;
    options_cur.optTol = optTol;
    options_cur.maxIter = maxIter;
    options_cur.verbose = verbose;
    options_cur.initType = initType;
    
    %[x, f] = TPower(A, options_cur);  
    [x, f] = TPower(A, options_cur);  
    
    X = [X,x];
    F = [F, f];    
    % deflation step
    if(m>1)
        P = eye(dim)- x*x';
        A = P*A*P;
    end
    
end

end


