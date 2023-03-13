 %% ------------------------------------------------------------------------
%%  Dynamic Prediction Pooling 
%% ------------------------------------------------------------------------
%    clear all;
%    clc;
% 
 tic
 varforecast = 6;
k = varforecast+1; 


     parameter = 'all'    %% estimate three parameters: rho, mu, sigma 
% parameter = 'one'   %% estimate one parameters: only rho

disp(' Start Dynamic Prediction Pool')

h=8
nperticle   = 500;
 nsim        = 5000;
 nburn       = ceil(0.3*nsim);
 
switch parameter
    case 'all'
      if h==8
        horizon = 8;
        load ./data/DP_WW_mean1_8.mat        % YY_mean YY_median; 
        load ./data/save_filter_DP1_8.mat    % save_filter  lambdasim        
        load('./results/MS_sample_02-Sep-2018_8.mat')           
                 WW_mean     = mean(save_Wt(nsim-nburn:end,:),1)'; 
%                 WW_mean     = median(save_Wt(nsim-nburn:end,:),1)';
      elseif h==4 
        horizon = 4;
        load ./data/DP_WW_mean1_4.mat        % YY_mean YY_median; 
        load ./data/save_filter_DP1.mat    % save_filter  lambdasim        
        load('./results/MS_sample_02-Sep-2018_4.mat')           
                 WW_mean     = mean(save_Wt(nsim-nburn:end,:),1)'; 
%                 WW_mean     = median(save_Wt(nsim-nburn:end,:),1)';
      elseif h==2 
        horizon = 2;
        load ./data/DP_WW_mean1_2.mat        % YY_mean YY_median; 
        load ./data/save_filter_DP1_2.mat    % save_filter  lambdasim        
        load('./results/MS_sample_02-Sep-2018_2.mat')           
                 WW_mean     = mean(save_Wt(nsim-nburn:end,:),1)';     
%                 WW_mean     = median(save_Wt(nsim-nburn:end,:),1)';
        end         
    case 'one'
        load ./data/DP_WW_mean2.mat         % YY_mean YY_median;
        load ./data/save_filter_DP2.mat     % save_filter  lambdasim
end   

 [Tobs]   = size(WW_mean,1); 
% T0         = 4;

 
rho_old     = 0.2;
mu_old      = 0; % icdf('norm',0.6,0,1);
sigma_old   =  1.0;
T0_Forecast = horizon+1;
h_Forecast  = 0;
 


% ------------------------------------------------------------------------
%  Graph of Filtered Posterior lambda (Dynamic Prediction Pools)
% ------------------------------------------------------------------------
%  Dynamic Prediction Pools  
a         =   0.95;       % (1-2*rate) percentage of Credibile interval of Parameters  
rate_a    =   (1-a)/2;
b         =   0.68;       % (1-2*rate) percentage of Credibile interval of Parameters  
rate_b    =   (1-b)/2;

% m       =   3;          %  Filtered Posterior lambda
% b       =   a/2/m;      %  Filtered Posterior lambda
% col1    =   0.75;   
% col2    =   0.65;
% col3    =   0.2;

sort_filter = zeros( nsim-nburn, Tobs );

%  Sort of sampling forecast
for k = 1:1:Tobs
    sort_filter(:,k) = sort(save_filter(nburn+1:end,k),1);
end  

%  Graph of Forecast 
 YY_mean = mean(save_filter(nburn+1:end,:),1)';
%  YY_median = median(save_filter(nburn+1:end,:),1)';
% YY_band = zeros(Tobs,2*m);
%   Sort_forecast_1 = sort_forecast(:,:);
YY_median = sort_filter(ceil((nsim-nburn)*0.5),:)'; 

% for j = 1:m 
%     YY_band(:,j)        = sort_filter(ceil((nsim-nburn)*(rate+b*(j-1))),:)'; 
%     YY_band(:,2*m-j+1)  = sort_filter(ceil((nsim-nburn)*(1-(rate+b*(j-1)))),:)';
% end

% YY_f_save = YY_band(:,1);
% for j = 1:m-1
%     YY_f_save = [ YY_f_save YY_band(:,j+1)-YY_band(:,j) ];
% end  
% YY_f_save = [ YY_f_save YY_median-YY_band(:,m) YY_band(:,m+1)-YY_median ];
% for j = m:-1:2
%     YY_f_save = [ YY_f_save YY_band(:,2*m+2-j)-YY_band(:,2*m+1-j) ];
% end

BC = csvread('./data/BC_Japan.csv',1,1);  

ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs)/4;  

h_D =figure('Position',[20,20,900,300],'File','Fig_7_dp_h',...
             'Name','Filtered Posterior lambda (Weight on Finanical Friction)','Color','w');
title_name = {'Filtered Posterior lambda (Weight on Finanical Friction model)'};
subplot(1,1,1)
hold on 
% hh = area(ti,YY_f_save(T0_Forecast:Tobs-h_Forecast,:)  ) ; 
% 
% for j = 1:m    
%     set(hh(j),'FaceColor',[col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3)])    
%     set(hh(2*(m+1)-j),'FaceColor',[col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3)])
% end
% set(hh(m+1),'FaceColor',[1-col1 1-col2 1-col3])
% set(hh(1),'FaceColor',[1 1 1])    
% set(hh,'LineStyle','none')      % Set all to same value

hh = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[0.5,0.65,1.0]) ;  

l1=plot(ti, WW_mean(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','--','Color','b', 'LineWidth',3.0);
% l2=plot(ti, YY_median(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','black', 'LineWidth',2.5);
l2=plot(ti, YY_mean(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','black', 'LineWidth',2.5);
l3=plot(ti, sort_filter(ceil((nsim-nburn)*rate_b),T0_Forecast:Tobs-h_Forecast)','LineStyle',':','Color','k', 'LineWidth',1.5); 
l4=plot(ti, sort_filter(ceil((nsim-nburn)*(1-rate_b) ),T0_Forecast:Tobs-h_Forecast)','LineStyle',':','Color','k', 'LineWidth',1.5); 
l5=plot(ti, 0.5*ones(Tobs-h_Forecast-T0_Forecast+1,1),'LineStyle','-','Color','red', 'LineWidth',2.5); 

% title( title_name(1),'FontSize',14 )  
ylim([0.0 1.0]);
xlim([1980 2000]);
legend([l2,l3,l1,hh ],{'Mean','68% band','MS Model','Recession'},'FontSize',12 )
hold off


 


