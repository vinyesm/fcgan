function H = gradient(Z,inputData,param)

switch param.f
    case 1 % prox
        H = Z - inputData.Y;
    case 2 % multitask
        H = inputData.X'*(inputData.X*Z - inputData.Y);
    case 3 % quadratic regression
        H = inputData.X'*diag(diag(inputData.X*Z*inputData.X')-inputData.Y)*inputData.X;
    case 4 % bilinear
        H = inputData.X1'*(inputData.X1*Z*inputData.X2 - inputData.Y)*inputData.X2';
end

end