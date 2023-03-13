
 horizon=8

% load('./data/MS_sample_12-jul-2015.mat')
 load('./data/MS_sample_19-May-2016.mat')
   MS_median= YY_median;
load(strcat('./data/save_filter_DP1_', num2str(horizon) , '.mat'))
% load('./data/save_filter_DP2.mat')
 YY_mean = YY_median;

nsim = 100;
nburn       = ceil(0.3*nsim);
T0_Forecast = horizon+1;
h_Forecast  = 0;

%  Dynamic Prediction Pools  
a       =   0.90;       % (1-2*rate) percentage of Credibile interval of Parameters  
rate    =   (1-a)/2;
m       =   3;          %  Filtered Posterior lambda
b       =   a/2/m;      %  Filtered Posterior lambda
col1    =   0.75;   
col2    =   0.65;
col3    =   0.2;
Tobs= 72;
sort_filter = zeros( nsim-nburn, Tobs );

%  Sort of sampling forecast
for k = 1:1:Tobs
    sort_filter(:,k) = sort(save_filter(nburn+1:end,k),1);
end  

%  Graph of Forecast 
 YY_mean = mean(save_filter(nburn+1:end,:),1)';
%  YY_median = median(save_filter(nburn+1:end,:),1)';
YY_band = zeros(Tobs,2*m);
%   Sort_forecast_1 = sort_forecast(:,:);
YY_median = sort_filter(ceil((nsim-nburn)*0.5),:)'; 

for j = 1:m 
    YY_band(:,j)        = sort_filter(ceil((nsim-nburn)*(rate+b*(j-1))),:)'; 
    YY_band(:,2*m-j+1)  = sort_filter(ceil((nsim-nburn)*(1-(rate+b*(j-1)))),:)';
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

h_D =figure('Position',[20,20,900,600],'Name','Filtered Posterior lambda (Weight on Finanical Friction)','Color','w');
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
set(hh,'LineStyle','none')        % Set all to same value
plot(ti, MS_median(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','--','Color','r', 'LineWidth',2.5);
plot(ti, YY_median(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','black', 'LineWidth',2.5); 
title( title_name(1),'FontSize',14 )  
hold off

% figure('Name','Predictive Densities of FF, DSGE','Color','w');   
% subplot(1,1,1)
%     plot(ti', log(PredDen_SW(T0_Forecast:Tobs-h_Forecast)), 'r', 'LineWidth',2)
% hold on
%     plot(ti', log(PredDen_KK(T0_Forecast:Tobs-h_Forecast)), 'b--', 'LineWidth',2) 
% %     plot(ti', log(PredDen_OP(T0_Forecast:Tobs-h_Forecast)),'g-', 'LineWidth',2) 
% hold off
%         
% title( 'Predictive Densities of New Keynesian and Finanical Friction','FontSize',12 );
% legend('New Keynesian', 'Financial Friction');
