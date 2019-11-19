
clc;clear;close all

aCd0 = [1;1e-4];
year = 2003;
opts = optimset('PlotFcn',@optimplotfval);
[x,ffinal] = fminsearch(@(aCd) optimize_alpha_CD_ref(aCd,year),aCd0,opts)
% non-smooth [alpha;cd]
% 1.270451301440938
%    0.001688020331327
% ffinal = 29.908

aCd0 = [1;1e-4];
year = 2004;
opts = optimset('PlotFcn',@optimplotfval);
[x,ffinal] = fminsearch(@(aCd) optimize_alpha_CD_ref(aCd,year),aCd0,opts)
% non-smooth [alpha;cd]


aCd0 = [1;1e-4];
year = 2005;
opts = optimset('PlotFcn',@optimplotfval);
[x,ffinal] = fminsearch(@(aCd) optimize_alpha_CD_ref(aCd,year),aCd0,opts)
% non-smooth [alpha;cd]
% 


aCd0 = [1;1e-4];
year = 2006;
opts = optimset('PlotFcn',@optimplotfval);
[x,ffinal] = fminsearch(@(aCd) optimize_alpha_CD_ref(aCd,year),aCd0,opts)
% non-smooth [alpha;cd]


% opts = optimset('PlotFcn',@optimplotfval);
% aCd0 = [1;1e-3];
% year = 2007;
% [x,ffinal] = fminsearch(@(aCd) optimize_alpha_CD_ref(aCd,year),aCd0,opts)
% x =
% 
%    2.663021548847048
%    0.001554770808786
% 
% 
% ffinal =
% 
%   27.296087370951255






