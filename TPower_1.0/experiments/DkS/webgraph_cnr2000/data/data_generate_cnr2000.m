clear all;

try
    G = bvgraph('cnr-2000');
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
W = ( W + W' ) / 2;
save ('data.mat', 'W');