function startingZ = set_default_Z(inputData,param)
if param.f == 1
    startingZ = zeros(size(inputData.Y));
elseif param.f == 2
    startingZ = zeros(size(inputData.X,2), size(inputData.Y,2));
elseif param.f == 3
    startingZ = zeros(size(inputData.X,2));
elseif param.f == 4
    startingZ = zeros(size(inputData.X1,2), size(inputData.X2,1));
end
end