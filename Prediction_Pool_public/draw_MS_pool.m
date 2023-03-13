% ------------------------------------------------------------------------
%   Markov Switching Predition Pooling 
% ------------------------------------------------------------------------
 clear all;
 clc;
 
 load('./data/MS_sample_17-Aug-2018_8.mat')

 horizon = 8;
 
varforecast = 6;
k = varforecast+1; 
disp(' Start Markov Switching Prediction Pool')

path('function',path);
 path('..\DSGE_KK2\results\FF2011',path);
    load(strcat('FF2011pred_den_', num2str(horizon), '.mat'))
       PredDen_KK = pred_den(:,k);
 path('..\DSGE_SW\results\sw2011',path);     
    load( strcat('sw2011pred_den_', num2str(horizon), '.mat'))
       PredDen_SW = pred_den(:,k);
%     load 'OP2011pred_den.mat'
%        PredDen_OP = pred_den;

% %pd = makedist('Normal');
[Tobs,n]    = size(PredDen_SW); 
T0_Forecast = horizon+1;
h_Forecast  = 0;
%   h_Forecast  = Tobs - 23; % pre_bubble
 
nsim    = 100000;
nburn   = round(0.4*nsim);

likeold     = -10^(1000);
priornew    = 0;
priorold    = 0;
St_old      = zeros(Tobs,1);
lambdasim   = zeros(nsim,2);

%*********************************************
%
%  Graph of Prob of Regime (Dynamic Prediction Pools)
%
%*********************************************

a       =   0.90;       % (1-2*rate) percentage of Credibile interval of Parameters  
rate    =   (1-a)/2;
m       =   3;          %  Filtered Posterior lambda
b       =   a/2/m;      %  Filtered Posterior lambda
col1    =   0.75;   
col2    =   0.65;
col3    =   0.2;

sort_filter = zeros( nsim-nburn, Tobs );
sort_weight = zeros( nsim-nburn, Tobs );

%  sort of sampling forecast
for k = 1:1:Tobs 
    sort_filter(:,k) = sort(save_St(nburn+1:end,k),1);  
    sort_weight(:,k) = sort(save_Wt(nburn+1:end,k),1);  
end

%  Graph of Forecast
YY_mean = mean(save_St(nsim-nburn:end,:),1)';
YY_band = zeros(Tobs,2*m);
%   sort_forecast_1 = sort_forecast(:,:);
YY_median = sort_filter(round((nsim-nburn)*0.5),:)'; 
%   YY_mean = mean(sort_filter(round((nsim-nburn)*(rate+b*(m+1))):round((nsim-nburn)*(1-(rate+b*(m+1)))),:),1)'; 

for j = 1:m 
    YY_band(:,j) = mean(sort_filter(1:round((nsim-nburn)*(1-(rate+b*(m-j+1)))),:),1)'; 
    YY_band(:,2*m-j+1) =mean( sort_filter(round((nsim-nburn)*(rate+b*(m-j+1))):nsim-nburn,:),1 )';
end
YY_f_save = YY_band(:,1);
for j = 1:m-1
    YY_f_save = [ YY_f_save YY_band(:,j+1)-YY_band(:,j) ];
end  
YY_f_save = [ YY_f_save YY_median-YY_band(:,m) YY_band(:,m+1)-YY_median ];
for j = m:-1:2
    YY_f_save = [ YY_f_save YY_band(:,2*m+2-j)-YY_band(:,2*m+1-j) ];
end

% ------------------------------------------------------------------------
%  Calculating Weight
% ------------------------------------------------------------------------
WW_mean     = mean(save_Wt(nsim-nburn:end,:),1)';  
WW_band     = zeros(Tobs,2*m);
WW_median   = sort_weight(round((nsim-nburn)*0.5),:)'; 
%   sort_weight_1 = squeeze(sort_weight(:,:,i)); 
 
for j = 1:m 
    WW_band(:,j) = sort_weight(round((nsim-nburn)*(rate+b*(j-1))),:)'; 
    WW_band(:,2*m-j+1) = sort_weight(round((nsim-nburn)*(1-(rate+b*(j-1)))),:)';
end
WW_save = WW_band(:,1);

for j = 1:m-1
    WW_save = [ WW_save WW_band(:,j+1)-WW_band(:,j) ];
end  
WW_save = [ WW_save WW_median-WW_band(:,m) WW_band(:,m+1)-WW_median ];
for j = m:-1:2
    WW_save = [ WW_save WW_band(:,2*m+2-j)-WW_band(:,2*m+1-j) ];
end
%    
% ------------------------------------------------------------------------  
%  drawing graph  
% ------------------------------------------------------------------------  

ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs-h_Forecast)/4; 

% h_MS2 = figure('Position',[20,20,900,600],'Name','Prob of Regime 1','Color','w');
% title_name={ '(a) Probabilities of Regime that FF model is better, \lambda_2','FontSize',14};
% subplot(2,1,1);
% hold on 
% hh = area(ti,  YY_f_save(T0_Forecast:Tobs-h_Forecast,:)) ; 
%       
% for j = 1:m    
%     set(hh(j),'FaceColor',...
%     [col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3) ])    
%     set(hh(2*(m+1)-j),'FaceColor',...
%     [col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3) ])
% end
% set(hh(m+1),'FaceColor',[1-col1 1-col2 1-col3 ])
% set(hh(1),'FaceColor',[1 1 1 ])    
% set(hh,'LineStyle','none') % Set all to same value
%       
% ylim([0.0 1.0]);  
%    plot(ti, YY_median(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','--','Color','r', 'LineWidth',2.5);
%    plot(ti, YY_mean(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','k', 'LineWidth',2.5);
% title( title_name(1),'FontSize',14 )  
% hold off
% 
% subplot(2,1,2);
% hold on 
% hh_w = area(ti, WW_save(T0_Forecast:Tobs-h_Forecast,:)) ; 
%       
% for j = 1:m    
%     set(hh_w(j),'FaceColor',...
%     [col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3) ])    
%     set(hh_w(2*(m+1)-j),'FaceColor',...
%     [col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3) ])
% end
% set(hh_w(m+1),'FaceColor',[1-col1 1-col2 1-col3 ])
% set(hh_w(1),'FaceColor',[1 1 1 ])    
% set(hh_w,'LineStyle','none')
% 
% %lambda_mean =  lambda0*(1-YY_mean(T0:Tobs-h,1))+lambda1*(YY_mean(T0:Tobs-h,1));
%           
%  plot(ti, WW_median(T0_Forecast:Tobs-h_Forecast,1), 'LineStyle','--','Color','r', 'LineWidth',2.5);   
%  plot(ti, WW_mean(T0_Forecast:Tobs-h_Forecast,1), 'LineStyle','-','Color','k', 'LineWidth',2.5); 
% ylim([0.0 1.0]); 
% title( '(b) Posterior Model Weight on FF Model, \lambda_t','FontSize',14 );
% hold off
% 
%  save \DSGE\DSGE_KK/Prediction_Pool/data/MS_WW_mean.mat  WW_mean WW_median  ; 
%  
 
%
 
BC = csvread('./data/BC_Japan.csv',1,1);  

h_MS3 = figure('Position',[20,20,900,600],'Name','Prob of Regime 1','Color','w');
title_name={ '(a) Probabilities of Regime that FF model is better; Prob(S_t = S_2)','FontSize',14};
subplot(2,1,1);
hold on 
% hh = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[1.0,0.75,1.0]) ;
hh = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[0.65,0.99,0.65]) ;
      
ylim([0.0 1.0]);  
   plot(ti, YY_median(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','b', 'LineWidth',2.5);
   plot(ti, YY_mean(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','--','Color','k', 'LineWidth',2.5);
   plot(ti, 0.5*ones(Tobs-h_Forecast-T0_Forecast+1,1),'LineStyle','-','Color','red', 'LineWidth',1.5); 
title( title_name(1),'FontSize',14 ) 
legend('Recession','Median','Mean')
hold off

subplot(2,1,2);
hold on 
% hh_w = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[1.0,0.75,1.0]) ;
hh_w = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[0.65,0.99,0.65]) ;
         
 plot(ti, WW_median(T0_Forecast:Tobs-h_Forecast,1), 'LineStyle','-','Color','b', 'LineWidth',2.5);   
 plot(ti, WW_mean(T0_Forecast:Tobs-h_Forecast,1), 'LineStyle','--','Color','k', 'LineWidth',2.5); 
 plot(ti, 0.5*ones(Tobs-h_Forecast-T0_Forecast+1,1),'LineStyle','-','Color','red', 'LineWidth',1.5); 
ylim([0.0 1.0]); 
title( '(b) Posterior Model Weight on FF Model; \lambda_t','FontSize',14 );
legend('Recession','Median','Mean')
hold off
 

