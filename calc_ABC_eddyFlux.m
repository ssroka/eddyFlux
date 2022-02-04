
salinity = 34*ones(m,n);% ppt for Lv calculation

for i = 1:length(year_vec)
    year = year_vec(i);
    dataFile = sprintf('%sERA5_patch_data_%d.mat',data_src,year);
    load(dataFile)
        
    % coefficients
    load(sprintf('opt_abCD_%sfilt_%s_L_%d_box%d_%d_%s_%d',con_str,filter_type,L/1000,box_num,er_box_num,model_str,year_vec(i)),'abCD','FFINAL','box_opt');
    filename = sprintf('Qs_QL_optimization_data_L_%d_filt_%s_box%d_%s_%d',L/1000,filter_type,box_num,model_str,year);

    p = size(SST_patch,3); % re calculate for leap year
    
    ZEROS = zeros(m,n,p);
    
    switch model_str
        case 'alpha'
            
            abCD = abCD.*abCD_factor;
            
            a = abCD(1);
            CD =  abCD(2);
            
            A = ZEROS;
            B = ZEROS;
            C1 = ZEROS;
            C2 = ZEROS;
            C3 = ZEROS;
            D = ZEROS;
            E1 = ZEROS;
            E2 = ZEROS;
            sshf_patch_bar = ZEROS;
            slhf_patch_bar = ZEROS;
            adh_CTRL = ZEROS;
            adh_prime = ZEROS;
            U_CTRL = ZEROS;
            U_prime = ZEROS;
            To_prime = ZEROS;
            dh_CTRL =  ZEROS;
            dh_prime =  ZEROS;
            qs_prime =  ZEROS;
            SST =  ZEROS;
            
            count = 1;
            fprintf('\n')
            for tt = 1:intvl:p % time points
                fprintf(' processing snapshot %d of %d\n',tt,p)
                
                %             debug_flag = false;
                
                % SST
                [SST_patch_CTRL,SST_prime] = FFT2D_filter(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                
                % qo
                [qo_CTRL,qo_prime] = FFT2D_filter(qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                %             SST_prime = SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt) - SST_patch_CTRL;
                
                Lv = SW_LatentHeat(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),'K',salinity,'ppt');
                
                % h
                h_diff = Lv.*(qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt)-qa_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt)) +...
                    c_p_air.*(DT_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt));
                [h_diff_CTRL,h_diff_prime] = FFT2D_filter(h_diff,dx,cf,debug_flag,lat_box,lon_box);
                %             h_diff_prime = h_diff - h_diff_CTRL;
                
                % alpha h
                ah_diff = a.*Lv.*(qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt)-qa_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt)) +...
                    a.*c_p_air.*(DT_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt));
                [ah_diff_CTRL,ah_diff_prime] = FFT2D_filter(ah_diff,dx,cf,debug_flag,lat_box,lon_box);
                %             ah_diff_prime = ah_diff - ah_diff_CTRL;
                
                % U
                [U_mag_CTRL,U_mag_prime] = FFT2D_filter(U_mag(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                %             U_mag_prime = U_mag(lon_patch_2_box_TF,lat_patch_2_box_TF,tt) - U_mag_CTRL;
                
                A_full = rho_a*CD*U_mag_CTRL.*h_diff_CTRL;
                B_full = rho_a*CD*U_mag_prime.*h_diff_prime;
                C1_full = rho_a*CD*U_mag_CTRL.*SST_prime.*ah_diff_prime;
                C2_full = rho_a*CD*U_mag_prime.*SST_prime.*ah_diff_CTRL;
                C3_full = rho_a*CD*U_mag_prime.*SST_prime.*ah_diff_prime;
                D_full = rho_a*CD*U_mag_prime.*h_diff_CTRL;
                E1_full = rho_a*CD*SST_prime.*U_mag_CTRL.*ah_diff_CTRL;
                E2_full = rho_a*CD*U_mag_CTRL.*h_diff_prime;
                
                
                if fft_first_flag
                    %                 debug_flag = 1;
                    A(:,:,count) = FFT2D_filter(A_full,dx,cf,debug_flag,lat_box,lon_box);
                    B(:,:,count) = FFT2D_filter(B_full,dx,cf,debug_flag,lat_box,lon_box);
                    C1(:,:,count) = FFT2D_filter(C1_full,dx,cf,debug_flag,lat_box,lon_box);
                    C2(:,:,count) = FFT2D_filter(C2_full,dx,cf,debug_flag,lat_box,lon_box);
                    C3(:,:,count) = FFT2D_filter(C3_full,dx,cf,debug_flag,lat_box,lon_box);
                    D(:,:,count) = FFT2D_filter(D_full,dx,cf,debug_flag,lat_box,lon_box);
                    E1(:,:,count) = FFT2D_filter(E1_full,dx,cf,debug_flag,lat_box,lon_box);
                    E2(:,:,count) = FFT2D_filter(E2_full,dx,cf,debug_flag,lat_box,lon_box);
                    sshf_patch_bar(:,:,count) = FFT2D_filter(sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                    slhf_patch_bar(:,:,count) = FFT2D_filter(slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                    fft_str = '';
                    
                else
                    A(:,:,count) = A_full;
                    B(:,:,count) = B_full;
                    C1(:,:,count) = C1_full;
                    C2(:,:,count) = C2_full;
                    C3(:,:,count) = C3_full;
                    D(:,:,count) = D_full;
                    E1(:,:,count) = E1_full;
                    E2(:,:,count) = E2_full;
                    sshf_patch_bar(:,:,count) = sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt);
                    slhf_patch_bar(:,:,count) = slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt);
                    
                    fft_str = 'noFFT_';
                end
                
                
                dh_CTRL(:,:,count) = h_diff_CTRL;
                dh_prime(:,:,count) = h_diff_prime;
                adh_CTRL(:,:,count) = ah_diff_CTRL;
                adh_prime(:,:,count) = ah_diff_prime;
                U_CTRL(:,:,count) = U_mag_CTRL;
                U_prime(:,:,count) = U_mag_prime;
                To_prime(:,:,count) = SST_prime;
                qs_prime(:,:,count) = qo_prime;
                SST(:,:,count) = SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt);
                
                
                count = count + 1;
            end
            
          

            save(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%d_%s_%d%s',L/1000,con_str,fft_str,filter_type,box_num,er_box_num,model_str,year,abCD_fac_str),...
                'A','B','C1','C2','C3','D','E1','E2','slhf_patch_bar','sshf_patch_bar','SST',...
                'dh_CTRL','dh_prime','adh_CTRL','adh_prime','U_CTRL','U_prime','To_prime','qs_prime',...
                'dx','cf','debug_flag','lat_box','lon_box','a','CD')
            fprintf('\nsaving A,B and C for %d\n',year)
            
            
        case 'beta'
            
            abCD = abCD.*abCD_factor;
            
            b = abCD(1);
            CD = abCD(2);
            
            
            Ab = ZEROS;
            Bb = ZEROS;
            Cb = ZEROS;
            Db = ZEROS;
            sshf_patch_bar = ZEROS;
            slhf_patch_bar = ZEROS;
            U_bar = ZEROS;
            U_prime = ZEROS;
            To_prime = ZEROS;
            dh_bar =  ZEROS;
            dh_prime =  ZEROS;
            Dq_prime =  ZEROS;
            SST =  ZEROS;
            
            load(filename,'bs_multiplier','bL_multiplier','U_bar','SST_prime','sshf_patch','slhf_patch')

            
            count = 1;
            fprintf('\n')
            for tt = 1:intvl:p % time points
                fprintf(' processing snapshot %d of %d\n',tt,p)
                
                %             debug_flag = false;
                
                % SST
                [SST_bar_tt,SST_prime_tt] = FFT2D_filter(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                
                % DT
                [DT_bar_tt,DT_prime_tt] = FFT2D_filter(DT_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                
                % qo
                [qo_CTRL_tt,qo_prime_tt] = FFT2D_filter(qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                
                % Dq
                Dq = (qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt)-qa_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt));
                [Dq_bar_tt,Dq_prime_tt] = FFT2D_filter(Dq,dx,cf,debug_flag,lat_box,lon_box);
                
                Lv = SW_LatentHeat(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),'K',salinity,'ppt');
                
                % h
                CD_h_diff = CD*(bs_multiplier(:,:,tt) + bL_multiplier(:,:,tt));
                %             CD_h_diff = rho_a*CD.*Lv.*Dq +...
                %                 rho_a*CD.*c_p_air.*(DT_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt));
                [CD_h_diff_bar_tt,CD_h_diff_prime_tt] = FFT2D_filter(CD_h_diff,dx,cf,debug_flag,lat_box,lon_box);
                %             h_diff_prime = h_diff - h_diff_CTRL;
                
                
                % U
                [U_bar_tt,U_prime_tt] = FFT2D_filter(U_mag(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                %             U_mag_prime = U_mag(lon_patch_2_box_TF,lat_patch_2_box_TF,tt) - U_mag_CTRL;
                
                Ab_full = U_bar_tt.*CD_h_diff;
                Bb_full = U_bar_tt.*CD_h_diff_prime_tt;
                Cb_full = b.*SST_prime_tt.*CD_h_diff_prime_tt;
                Db_full = b.*SST_prime_tt.*CD_h_diff_bar_tt;
                
                
                if fft_first_flag
                    %                 debug_flag = 1;
                    Ab(:,:,count) = FFT2D_filter(Ab_full,dx,cf,debug_flag,lat_box,lon_box);
                    Bb(:,:,count) = FFT2D_filter(Bb_full,dx,cf,debug_flag,lat_box,lon_box);
                    Cb(:,:,count) = FFT2D_filter(Cb_full,dx,cf,debug_flag,lat_box,lon_box);
                    Db(:,:,count) = FFT2D_filter(Db_full,dx,cf,debug_flag,lat_box,lon_box);
                    sshf_patch_bar(:,:,count) = FFT2D_filter(sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                    slhf_patch_bar(:,:,count) = FFT2D_filter(slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                    fft_str = '';
                    
                else
                    A(:,:,count) = Ab_full;
                    Bb(:,:,count) = Bb_full;
                    Cb(:,:,count) = Cb_full;
                    Db(:,:,count) = Db_full;
                    
                    sshf_patch_bar(:,:,count) = sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt);
                    slhf_patch_bar(:,:,count) = slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt);
                    
                    fft_str = 'noFFT_';
                end
                
                
                dh_bar(:,:,count) = CD_h_diff_bar_tt./CD;
                dh_prime(:,:,count) = CD_h_diff_prime_tt./CD;
                U_bar(:,:,count) = U_bar_tt;
                U_prime(:,:,count) = U_prime_tt;
                To_prime(:,:,count) = SST_prime_tt;
                Dq_prime(:,:,count) = Dq_prime_tt;
                SST(:,:,count) = SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt);
                
                
                count = count + 1;
            end
            save(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%d_%s_%d%s',L/1000,con_str,fft_str,filter_type,box_num,er_box_num,model_str,year,abCD_fac_str),...
                'Ab','Bb','Cb','Db','slhf_patch_bar','sshf_patch_bar','SST',...
                'dh_bar','dh_prime','U_bar','U_prime','To_prime','Dq_prime',...
                'dx','cf','debug_flag','lat_box','lon_box','b','CD')
            fprintf('\nsaving Ab,Bb and Cb for %d\n',year)
            
            
        case 'alphabeta'
            
            abCD = abCD.*abCD_factor;
            
            a = abCD(1);
            b = abCD(2);
            CD = abCD(3);
            
            Aab = ZEROS;
            Bab = ZEROS;
            C1ab = ZEROS;
            C2ab = ZEROS;
            D1ab = ZEROS;
            D2ab = ZEROS;
            E1ab = ZEROS;
            E2ab = ZEROS;
            sshf_patch_bar = ZEROS;
            slhf_patch_bar = ZEROS;
            U_bar = ZEROS;
            U_prime = ZEROS;
            To_prime = ZEROS;
            dh_bar =  ZEROS;
            dh_prime =  ZEROS;
            Dq_prime =  ZEROS;
            SST =  ZEROS;
            
            load(filename,'abs_multiplier','abL_multiplier','U_bar','SST_prime','sshf_patch','slhf_patch')

            
            count = 1;
            fprintf('\n')
            for tt = 1:intvl:p % time points
                fprintf(' processing snapshot %d of %d\n',tt,p)
                
                %             debug_flag = false;
                
                % SST
                [SST_bar_tt,SST_prime_tt] = FFT2D_filter(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                
                % DT
                [DT_bar_tt,DT_prime_tt] = FFT2D_filter(DT_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                
                % qo
                [qo_CTRL_tt,qo_prime_tt] = FFT2D_filter(qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                
                % Dq
                Dq = (qo_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt)-qa_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt));
                [Dq_bar_tt,Dq_prime_tt] = FFT2D_filter(Dq,dx,cf,debug_flag,lat_box,lon_box);
                
                Lv = SW_LatentHeat(SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),'K',salinity,'ppt');
                
                % h
                CD_h_diff = CD*(abs_multiplier(:,:,tt) + abL_multiplier(:,:,tt));
                %             CD_h_diff = rho_a*CD.*Lv.*Dq +...
                %                 rho_a*CD.*c_p_air.*(DT_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt));
                [CD_h_diff_bar_tt,CD_h_diff_prime_tt] = FFT2D_filter(CD_h_diff,dx,cf,debug_flag,lat_box,lon_box);
                %             h_diff_prime = h_diff - h_diff_CTRL;
                
                
                % U
                [U_bar_tt,U_prime_tt] = FFT2D_filter(U_mag(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                %             U_mag_prime = U_mag(lon_patch_2_box_TF,lat_patch_2_box_TF,tt) - U_mag_CTRL;
                
                Aab_full = U_bar_tt.*CD_h_diff;
                Bab_full = U_bar_tt.*CD_h_diff_prime_tt;
                C1ab_full = a.*SST_prime_tt.*U_bar_tt.*CD_h_diff;
                C2ab_full = a.*SST_prime_tt.*U_bar_tt.*CD_h_diff_prime_tt;
                D1ab_full = b.*SST_prime_tt.*CD_h_diff;
                D2ab_full = b.*SST_prime_tt.*CD_h_diff_prime_tt;
                E1ab_full = a.*b.*SST_prime_tt.*SST_prime_tt.*CD_h_diff;
                E2ab_full = a.*b.*SST_prime_tt.*SST_prime_tt.*CD_h_diff_prime_tt;
                
                
                if fft_first_flag
                    %                 debug_flag = 1;
                    Aab(:,:,count) = FFT2D_filter(Aab_full,dx,cf,debug_flag,lat_box,lon_box);
                    Bab(:,:,count) = FFT2D_filter(Bab_full,dx,cf,debug_flag,lat_box,lon_box);
                    C1ab(:,:,count) = FFT2D_filter(C1ab_full,dx,cf,debug_flag,lat_box,lon_box);
                    C2ab(:,:,count) = FFT2D_filter(C2ab_full,dx,cf,debug_flag,lat_box,lon_box);
                    D1ab(:,:,count) = FFT2D_filter(D1ab_full,dx,cf,debug_flag,lat_box,lon_box);
                    D2ab(:,:,count) = FFT2D_filter(D2ab_full,dx,cf,debug_flag,lat_box,lon_box);
                    E1ab(:,:,count) = FFT2D_filter(E1ab_full,dx,cf,debug_flag,lat_box,lon_box);
                    E2ab(:,:,count) = FFT2D_filter(E2ab_full,dx,cf,debug_flag,lat_box,lon_box);
                    sshf_patch_bar(:,:,count) = FFT2D_filter(sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                    slhf_patch_bar(:,:,count) = FFT2D_filter(slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt),dx,cf,debug_flag,lat_box,lon_box);
                    fft_str = '';
                    
                else
                    Aab(:,:,count) = Aab_full;
                    Bab(:,:,count) = Bab_full;
                    C1ab(:,:,count) = C1ab_full;
                    C2ab(:,:,count) = C2ab_full;
                    D1ab(:,:,count) = D1ab_full;
                    D2ab(:,:,count) = D2ab_full;
                    E1ab(:,:,count) = E1ab_full;
                    E2ab(:,:,count) = E2ab_full;
                    sshf_patch_bar(:,:,count) = sshf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt);
                    slhf_patch_bar(:,:,count) = slhf_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt);
                    
                    fft_str = 'noFFT_';
                end
                
                
                dh_bar(:,:,count) = CD_h_diff_bar_tt./CD;
                dh_prime(:,:,count) = CD_h_diff_prime_tt./CD;
                U_bar(:,:,count) = U_bar_tt;
                U_prime(:,:,count) = U_prime_tt;
                To_prime(:,:,count) = SST_prime_tt;
                Dq_prime(:,:,count) = Dq_prime_tt;
                SST(:,:,count) = SST_patch(lon_patch_2_box_TF,lat_patch_2_box_TF,tt);
                
                
                count = count + 1;
            end
            save(sprintf('ABC_terms_%d_%sfilt_%s%s_box%d_%d_%s_%d%s',L/1000,con_str,fft_str,filter_type,box_num,er_box_num,model_str,year,abCD_fac_str),...
                'Aab','Bab','C1ab','C2ab','D1ab','D2ab','E1ab','E2ab',...
                'slhf_patch_bar','sshf_patch_bar','SST',...
                'dh_bar','dh_prime','U_bar','U_prime','To_prime','Dq_prime',...
                'dx','cf','debug_flag','lat_box','lon_box','a','b','CD')% just added the a, not save yet
            fprintf('\nsaving Ab,Bb and Cb for %d\n',year)
            
            
            
            
            
            
    end
    
    
    
    
    
end
