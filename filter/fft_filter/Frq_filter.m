% =========================================================================
% Frq_filter
% =========================================================================
% Created By    : Sydney Sroka
% Last Edited By: Sydney Sroka
% Date          : 03/29/2021
%
% ------------------------
% This script was written to practice filtering a 1D signal

clear;close all;clc

% if you change the period length to anything that isn't a multiply of pi
% (or you could change the length of the input NT to some
% non-whole number) the output is not exact

T1 = 4; % number of cycles per second
T2 = 2; % number of cycles per second

% frq
f1 = T1/(2*pi);
f2 = T2/(2*pi);

N = 500; % length of signal

NT = 100; % number of seconds

t_vec = linspace(0,1,N);

cT = 1.5*pi; % cutoff period;

%

s = cos(f1*t_vec);
plot(t_vec,s,'linewidth',2)
