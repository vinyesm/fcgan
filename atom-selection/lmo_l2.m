function [maxval,new_atom]=lmo_l2(s,param)
maxval=norm('fro');
new_atom=s/maxval;
