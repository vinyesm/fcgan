function param = set_default_param(param)
%% set default values to param 



if ~isfield(param,'q')
    param.q = param.k;
end


if ~isfield(param,'PSD')
    param.PSD = false;
end

if ~isfield(param,'nbMainLoop')
    param.nbMainLoop = 100;
end


if ~isfield(param,'powerIter')
    param.powerIter = 100;
end
if ~isfield(param,'stPtPowerIter')
    param.stPtPowerIter = 100;
end

if ~isfield(param,'epsStop')
    param.epsStop = .1;
end

if ~isfield(param,'innerLoopIter')
    param.innerLoopIter = 100;
end

if ~isfield(param,'niterPS')
    param.niterPS =  200;
end

if ~isfield(param,'PSdualityEpsilon')
    param.PSdualityEpsilon = 1e-3;
end

if ~isfield(param,'verbose')
    param.verbose = 0;
end

if ~isfield(param,'debug')
    param.debug = 0;
end

if ~isfield(param,'sloppy')
    param.sloppy = 0;
end

end