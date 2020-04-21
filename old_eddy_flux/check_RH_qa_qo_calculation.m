ccc

addpath('~/Documents/MATLAB/util/')


for year = 2003
    
    load(sprintf('/Users/ssroka/MIT/Research/eddyFlux/ERA5_data/ERA5_patch_data_%d',year),...
        'qo_patch','qa_patch','RH_patch','P0_patch','SST_patch','t2m_patch','d2m_patch')    
    
    % ---- qo ----
    qo_patch_check = SAM_qsatWater(SST_patch(:,:,end), P0_patch(:,:,end)) ;
    qo = qo_patch(:,:,end);
    qo_diff = nanmean(qo(:) - qo_patch_check(:));
    
    % ---- qa ----
    e_sat = SAM_psatWater(t2m_patch(:,:,end)); %Pa
    e = RH_patch(:,:,end)./100.*e_sat;
    r = 0.622 * e ./ (P0_patch(:,:,end)-e);
    qa_patch_check = r./(1+r);
    qa = qa_patch(:,:,end);
    qa_diff = nanmean(qa(:) - qa_patch_check(:));
    
    % ---- RH ----
    RH_diff = zeros(size(d2m_patch,3),1);
    e_diff = zeros(size(d2m_patch,3),1);
    for tt = 1:size(d2m_patch,3)
        RH_tt = RH_patch(:,:,tt);
        T = t2m_patch(:,:,tt);
        TD = d2m_patch(:,:,tt);
        RH_patch_check = 100.*(exp((17.625.*TD)./(243.04+TD))./exp((17.625.*T)./(243.04+T)));
        RH_diff(tt) = nanmean(RH_tt(:) - RH_patch_check(:));
        
        TC = T-273.15;
        es = 0.6112*exp((17.67.*TC)./(243.5+TC))*1000; % Pa
        e_sat = SAM_psatWater(t2m_patch(:,:,tt)); % Pa
        e_diff(tt) = nanmean(es(:)-e_sat(:));
    end
    
end
qo_diff;
qa_diff;
RH_diff;

figure
plot(e_diff)
return


TC = T-273.15;
TDC = TD-273.15;
es = 0.6112*exp((17.67.*TC)./(243.5+TC)); %kPa
% Bolton, D. 1980. The computation of equivalent potential temperature.
% Mon. Wea. Rev.. 108. 1046?1053.
e = 0.6112*exp((17.67.*TDC)./(243.5+TDC)); %kPa
c = 243.04;
b = 17.625;
RH3 = 100*exp((c*b*(TDC-TC))./(c+TC)./(c+TDC));

RH_2 = e./es*100;
figure
imagesc(RH_2-RH_tt)
colorbar
figure
imagesc(RH_patch_check-RH_tt)
colorbar
figure
subplot(2,4,1)
imagesc(RH_tt)
colorbar
subplot(2,4,2)
imagesc(RH_2)
colorbar
subplot(2,4,3)
imagesc(RH_patch_check)
colorbar
subplot(2,4,4)
imagesc(RH3)
colorbar
subplot(2,4,5)
imagesc(RH_2-RH_tt)
colorbar
subplot(2,4,6)
imagesc(RH3-RH_tt)
colorbar
subplot(2,4,7)
imagesc(RH_patch_check-RH_tt)
colorbar
subplot(2,4,8)
imagesc(RH3-RH_2)
colorbar


