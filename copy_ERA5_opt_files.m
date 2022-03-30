

for year =  2004:2018
    f1 = sprintf('opt_abCD_filt_fft_L_250_box3_beta_%d.mat',year);
    f2 = sprintf('opt_abCD_filt_fft_L_250_box3_3_beta_ERA5_%d.mat',year);
    copyfile(f1,f2)
end
%%
for abCD_fac_vec = [-10 -8 -5 -4 -2 -0.5 1 0.5 2 4 5 8 10];

if strcmp(model_str,'alphabeta')
    abCD_factor = [abCD_fac_vec abCD_fac_vec 1]';
else
    abCD_factor = [abCD_fac_vec 1]';
end

if all(abs(abCD_factor-1)<1e-10)
    abCD_fac_str = '';
else
    abCD_fac_str = ['_'];
    for abCD_ii = 1:length(abCD_factor)
        abCD_fac_str = [abCD_fac_str strrep(num2str(abCD_factor(abCD_ii)),'.','_') '_'];
    end
    abCD_fac_str = abCD_fac_str(1:end-1);
end

for year =  2004:2018
    f1 = sprintf('mean_rms_err_%s_%d%s.mat',model_str,year,abCD_fac_str);
    f2 = sprintf('mean_rms_err_%s_%s_%d%s.mat',model_str,reanalysis_src,year,abCD_fac_str);
    copyfile(f1,f2)
end
end