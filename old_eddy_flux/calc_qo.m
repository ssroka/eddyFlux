function [qo] = calc_qo(sst,p,DT)

if nargin<2
    DT = 0.5; %K
end
addpath('/Users/ssroka/Documents/MATLAB/util/')
qo = SAM_qsatWater(sst-DT, p) ; 



end