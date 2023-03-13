 %% ------------------------------------------------------------------------
%%  Dynamic Prediction Pooling 
%% ------------------------------------------------------------------------
%    clear all;
%    clc;
% 
 tic
 varforecast = 6;
k = varforecast+1; 

%     parameter = 'all'    %% estimate three parameters: rho, mu, sigma 
%  parameter = 'one'   %% estimate one parameters: only rho

% disp(' Start Dynamic Prediction Pool')
% 
% path('function',path);
% path('c:\DSGE\DSGE_KK\DSGE_KK2\results\FF2011',path);
%  display('Loading new prediction density');
%     load 'FF2011pred_den.mat'
%        PredDen_KK = pred_den(:,k);
%  path('c:\DSGE\DSGE_KK\DSGE_SW\results\sw2011',path);      
%     load 'sw2011pred_den.mat'
%        PredDen_SW = pred_den(:,k);
% %     load 'OP2011pred_den.mat'
% %        PredDen_OP = pred_den;
% 
% switch parameter
%     case 'all'
%         load ./data/DP_WW_mean1.mat        % YY_mean YY_median; 
%         load ./data/save_filter_DP1.mat    % save_filter  lambdasim
%     case 'one'
%         load ./data/DP_WW_mean2.mat         % YY_mean YY_median;
%         load ./data/save_filter_DP2.mat     % save_filter  lambdasim
% end   

% ------------------------------------------------------------------------
%  Graph of Filtered Posterior lambda (Dynamic Prediction Pools)
% ------------------------------------------------------------------------


%  Sort of sampling forecast
% for k = 1:1:Tobs
%     sort_filter(:,k) = sort(save_filter(nburn+1:end,k),1);
% end  

%  Graph of Forecast 
%  YY_mean = mean(save_filter(nburn+1:end,:),1)';

BC = csvread('./data/BC_Japan.csv',1,1);  
rates = csvread('./data/loan_rates.csv',1,1);  

T0_Forecast = 5;
h_Forecast  = 0;

ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs)/4;  

h_D =figure('Position',[20,20,900,300],'Name','Filtered Posterior lambda (Weight on Finanical Friction)','Color','w');
title_name = {'Filtered Posterior lambda (Weight on Finanical Friction model)'};
subplot(1,1,1)
hold on 

hh = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[0.5,0.85,1.0]) ;  

plot(ti, rates(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','black', 'LineWidth',2.5); 
plot(ti, rates(T0_Forecast:Tobs-h_Forecast,2),'LineStyle',':','Color','b', 'LineWidth',2.5); 
plot(ti, rates(T0_Forecast:Tobs-h_Forecast,3),'LineStyle','--','Color','r', 'LineWidth',2.5); 

% title( title_name(1),'FontSize',14 )  
ylim([0.0 10]);
legend({'Recessions','policy rate','loan rate 1', 'loan rate 2'},'FontSize',12)
hold off

h_e =figure('Position',[20,20,900,300],'Name','Filtered Posterior lambda (Weight on Finanical Friction)','Color','w');
title_name = {'Filtered Posterior lambda (Weight on Finanical Friction model)'};
subplot(1,1,1)
hold on 

hh1 = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[0.5,0.85,1.0]) ;  
hh2 = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(-100),'LineStyle','non','FaceColor',[0.5,0.85,1.0]) ;

l2=plot(ti, rates(T0_Forecast:Tobs-h_Forecast,2)-rates(T0_Forecast:Tobs-h_Forecast,1), ...
    'LineStyle',':','Color','b', 'LineWidth',2.5); 
l3=plot(ti, rates(T0_Forecast:Tobs-h_Forecast,3)-rates(T0_Forecast:Tobs-h_Forecast,1),...
    'LineStyle','--','Color','r', 'LineWidth',2.5); 
plot(ti, 0*ones(Tobs-h_Forecast-T0_Forecast+1,1),'LineStyle','-','Color','k', 'LineWidth',1.5); 

% title( title_name(1),'FontSize',14 )  
ylim([-1 4]);
legend([hh1, l2, l3],{'Recessions', 'loan rate 1', 'loan rate 2'},'FontSize',12)
hold off


 


