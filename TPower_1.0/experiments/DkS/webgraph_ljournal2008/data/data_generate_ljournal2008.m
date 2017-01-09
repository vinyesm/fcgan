clear all;

try
    G = bvgraph('ljournal-2008');
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

save ('data.mat', 'W');

