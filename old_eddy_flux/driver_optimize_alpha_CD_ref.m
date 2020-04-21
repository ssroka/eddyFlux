clc;clear;close all
years = 2003:2007;
CD_ref_vec = [0.001688 0.00169 0.00158 0.00163 0.00155 0.001555];
alpha_vec =  [1.27045 1.5335 1.94 2.111 2.6630]*0.5;


options = optimoptions(@fmincon,'Display','iter','steptolerance',1e-5);
for i = 1:length(years)


% aCd0 = X{years(i)-2002};

% 6 params
% load X_FFINAL_2007.mat
% aCd0 = X{i};
% LB = [-Inf; 0; -Inf; 0; 0; 0];
% UB = [Inf; Inf; Inf; Inf; Inf; Inf];
% [x,ffinal] = fmincon(@(aCd) optimize_alpha_CD_ref(aCd,years(i)),aCd0,[],[],[],[],LB,UB,[],options);

% 5 params
% load X_FFINAL_2003_1CD.mat
% aCd0 = X{i};
% LB = [-Inf; 0; -Inf; 0; 0];
% UB = [Inf; Inf; Inf; Inf; Inf];
% [x,ffinal] = fmincon(@(aCd) optimize_alpha_CD_ref(aCd,years(i)),aCd0,[],[],[],[],LB,UB,[],options);


% 4 params
% load X_FFINAL_2007.mat
% aCd0 = aCd0([1 3 5 6]);
load(sprintf('X_FFINAL_%d_const_only',years(i)),'X','FFINAL');
aCd0 = X{i};

LB = [-Inf; -Inf; 0; 0];
UB = [Inf; Inf; Inf; Inf];
[x,ffinal] = fmincon(@(aCd) optimize_alpha_CD_ref(aCd,years(i)),aCd0,[],[],[],[],LB,UB,[],options);

X{i} = x;
FFINAL{i}  =ffinal;
save(sprintf('X_FFINAL_%d_const_only',years(i)),'X','FFINAL');
end








