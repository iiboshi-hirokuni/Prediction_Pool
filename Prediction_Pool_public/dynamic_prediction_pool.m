% ------------------------------------------------------------------------
%  Dynamic Prediction Pooling 
% ------------------------------------------------------------------------
clear all;
% clc;

horizon = 8;

varforecast = 6;
k = varforecast+1; 
%% 並列処理の　オン・オフ
% parallel = 'off' 
parallel = 'on'

disp(' Start Dynamic Prediction Pool')

path('function',path);
path('..\DSGE_KK2\results\FF2011',path);
 display('Loading new prediction density');
    load 'FF2011pred_den_8.mat'
       PredDen_KK = pred_den(:,k);
 path('..\DSGE\DSGE_KK\DSGE_SW\results\sw2011',path);      
    load 'sw2011pred_den_8.mat'
       PredDen_SW = pred_den(:,k);
%     load 'OP2011pred_den.mat'
%        PredDen_OP = pred_den;

[Tobs,n]   = size(PredDen_KK); 
T0         = 4;
nsim       = 5000;
rho        = 0.7;
mu         =  icdf('norm',0.5,0,1);
sigma      = 2.0;
T0_Forecast = horizon+1;
h_Forecast  = 0;

ncores = 4;
rv = randn(Tobs,1,nsim );

switch parallel
  case 'on'    
%      poolobj = parpool;
    tic
     [re_xt, log_lik] = ParticleFilter_par(PredDen_SW, PredDen_KK,...
                          T0_Forecast, h_Forecast, nsim, rho,mu,sigma, rv,ncores);
    toc
    disp(log_lik)
%      delete(poolobj)
  case 'off'
     tic
      [re_xt, log_lik] = ParticleFilter(PredDen_SW, PredDen_KK,...
                          T0_Forecast, h_Forecast, nsim, rho,mu,sigma);   
      disp(log_lik)
     toc
end

% ------------------------------------------------------------------------
%  Graph of Filtered Posterior lambda (Dynamic Prediction Pools)
% ------------------------------------------------------------------------
%  Dynamic Prediction Pools  
a       =   0.90;       % (1-2*rate) percentage of Credibile interval of Parameters  
rate    =   (1-a)/2;
m       =   3;          %  Filtered Posterior lambda
b       =   a/2/m;      %  Filtered Posterior lambda
col1    =   0.75;   
col2    =   0.65;
col3    =   0.2;

sort_filter = zeros( nsim, Tobs );

%  Sort of sampling forecast
for k = 1:1:Tobs
    sort_filter(:,k) = sort(re_xt(k,:)',1);
end  

%  Graph of Forecast 
YY_mean = mean(re_xt(:,:)',1)';
YY_band = zeros(Tobs,2*m);
%   Sort_forecast_1 = sort_forecast(:,:);
YY_median = sort_filter(round((nsim)*0.5),:)'; 

for j = 1:m 
    YY_band(:,j)        = sort_filter(round((nsim)*(rate+b*(j-1))),:)'; 
    YY_band(:,2*m-j+1)  = sort_filter(round((nsim)*(1-(rate+b*(j-1)))),:)';
end

YY_f_save = YY_band(:,1);
for j = 1:m-1
    YY_f_save = [ YY_f_save YY_band(:,j+1)-YY_band(:,j) ];
end  
YY_f_save = [ YY_f_save YY_median-YY_band(:,m) YY_band(:,m+1)-YY_median ];
for j = m:-1:2
    YY_f_save = [ YY_f_save YY_band(:,2*m+2-j)-YY_band(:,2*m+1-j) ];
end



ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs)/4;  

h_D = figure('Position',[20,20,900,600],'Name','Filtered Posterior lambda (Weight on Finanical Friction)','Color','w');
title_name = {'Filtered Posterior lambda (Weight on Finanical Friction model)'};
subplot(1,1,1)
hold on 
hh = area(ti,YY_f_save(T0_Forecast:Tobs-h_Forecast,:)  ) ; 

for j = 1:m    
    set(hh(j),'FaceColor',[col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3)])    
    set(hh(2*(m+1)-j),'FaceColor',[col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3)])
end
set(hh(m+1),'FaceColor',[1-col1 1-col2 1-col3])
set(hh(1),'FaceColor',[1 1 1])    
set(hh,'LineStyle','none')                                                      % Set all to same value
plot(ti, YY_mean(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','black', 'LineWidth',2.5); 
title( title_name(1),'FontSize',14 )  
hold off

figure('Name','Predictive Densities of FF, DSGE and Small Open','Color','w');   
subplot(1,1,1)
    plot(ti', log(PredDen_SW(T0_Forecast:Tobs-h_Forecast)), 'r', 'LineWidth',2)
hold on
    plot(ti', log(PredDen_KK(T0_Forecast:Tobs-h_Forecast)), 'b--', 'LineWidth',2) 
%     plot(ti', log(PredDen_OP(T0_Forecast:Tobs-h_Forecast)),'g-', 'LineWidth',2) 
hold off
        
title( 'Predictive Densities of New Keynesian, Finanical Friction and Small Open','FontSize',12 );
legend('New Keynesian', 'Financial Friction');

 save \DSGE\DSGE_KK/Prediction_Pool/data/DP_WW_mean.mat  YY_mean; 
 
 est_date = datestr(date);   
name = ['./results/Dynamic_Pool_',num2str(nsim),'_',est_date];
        saveas(h_D,name,'fig')
 
 
 