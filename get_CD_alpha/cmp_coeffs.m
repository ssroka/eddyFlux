clear;close all;clc


year_vec = [2003:2019];

L_vec = [500000];

filter_type = 'lanczos'; % filter type 'lanczos' or 'boxcar'

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
            
            load(sprintf('opt_aCD_%sfilt_%s_L_%d_%d_to_%d',con_str,filter_type,L/1000,year_vec(1),year_vec(y)),'X','FFINAL');

            plot(1:length(X{1}),X{year-2002},'color',c(y,:),'marker',m(i),'markerfacecolor',mrkr_clr,...
                'markersize',20,'displayname',sprintf('L = %d km%s, %d',L/1000,con_name,year),'linewidth',2)
            set(gca,'fontsize',20,'xtick',[1:4],'xticklabel',{'$$\alpha_s$$','$$\alpha_L$$','$$C_D^s$$','$$C_D^L$$'},'ticklabelinterpreter','latex')
            lh = legend('-dynamiclegend');
            set(lh,'interpreter','latex')
            hold on
            
        end
        
    end
    
    set(gcf,'color','w','position',[ 232  1  1188  801],'NumberTitle','off','Name',filter_type)
    
    

end

title(filter_type)
    update_figure_paper_size()

    print(sprintf('../imgs/cmp_alpha_CD_L_%d_%s',L/1000,filter_type),'-dpdf')







