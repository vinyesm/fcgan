clear all;

try
    G = bvgraph('hollywood-2009');
    fprintf('The library was already compiled!\n');
catch
    s = lasterror;
    if strcmp(s.identifier, 'bvgfun:notCompiled')
        fprintf('The library is now compiled!\n');
    else
        error('bvgraph:compileError','the library could not automatically compile!');
    end
end

W = sparse(G);

%W = ( W + W' ) / 2;
dim = size(W,1);
pieces = 100;
interval = floor(dim / pieces);

idx = [];
lb = 1;
ub = interval;
i = 1;
while (ub <= dim)

    idx_cur = [lb : ub];

    
    W_cur = W(idx_cur, :);
    idx = [idx, {idx_cur}];
   
    save (['data_piece_',int2str(i),'.mat'], 'W_cur');
    
    i = i+1;
    lb = (i-1)*interval + 1;
    ub = i*interval;
end

if (lb < dim)
    idx_cur = [lb : dim];
    W_cur = W(idx_cur, :);
    idx = [idx, {idx_cur}];
   
    save (['data_piece_',int2str(i),'.mat'], 'W_cur');
end

 save (['index.mat'], 'idx', 'dim');

