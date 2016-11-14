function [x,d,gamma]=as_tr_l1()

% [x,d,gamma] solve subproblem with bcmm
% s2 = gamma
% test <s2,ai><= mu forall i in Jc
% J <= J U {imax}, adding most violated constraint