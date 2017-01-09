function get_src

curdir = pwd;
try
    [filepath filename] = fileparts(mfilename('fullpath'));
    cd (filepath);
    
    type readme.txt
    
    input(['Press any key to download the GPL protected codes\n' ...
           'or Ctrl-C to stop here\n']);
    
    try
        urlwrite('http://wilkinson.stanford.edu/~dgleich/codes/bvgraph-1.2-src.zip',...
            'bvgraph-1.2-src.zip');
        
    catch
        fprintf(['Eek!  Getting the source failed!!\n' ...
            'Please send mithandor@gmail.com an email immediately\n']);
        fprintf('I''m so sorry, but is not going to compile without the source.\n');
        rethrow(lasterror);
    end
    
    unzip('bvgraph-1.2-src.zip');
    
    cd(curdir);
catch
    % make sure they end up back where they started
    cd(curdir);
    rethrow(lasterror);
end    
