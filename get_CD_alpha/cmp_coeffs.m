clear;close all;clc


year_vec = [2003:2019];

L_vec = [250000];

filter_type = 'fft'; % filter type 'lanczos' or 'boxcar'
param_num_str = '_3param';

m = 'so^';
c = [252 154 69;
    150 70 255;
    63 222 198]./256;

c=cool(length(year_vec));

for y = 1:length(year_vec)
    year = year_vec(y);
    for j = 1
        
        for i = 1:length(L_vec)
            
            if j==1
                con_name = '';
                con_str = '';
                mrkr_clr = 'none';
            else
                con_name = ', $$\alpha>0$$';
                con_str = 'cons_';
                mrkr_clr = c(i,:);
            end
            
            L = L_vec(i);
            
            load(sprintf('opt_aCD_%sfilt_%s_L_%d%s_%d_to_%d',con_str,filter_type,L/1000,param_num_str,year_vec(1),year_vec(y)),'X','FFINAL','box_opt');
            as(y) = X{year-2002}(1);
            aL(y) = X{year-2002}(2);
            CD(y) = X{year-2002}(3);
        end
    end
end

plot(year_vec,as,'kx-','displayname','$\alpha_s$','linewidth',2);
hold on
plot(year_vec,aL,'bs-','displayname','$\alpha_L$','linewidth',2);
xlabel('year')
lh = legend('-dynamiclegend');
set(lh,'interpreter','latex')
hold on
set(gca,'fontsize',25)
set(gcf,'color','w','position',[ 232  1  1188  801],'NumberTitle','off','Name',filter_type)

title('$\alpha$ coefficients','interpreter','latex')
update_figure_paper_size()

print(sprintf('../imgs/cmp_alpha_CD_L_%d_%s',L/1000,filter_type),'-dpdf')







