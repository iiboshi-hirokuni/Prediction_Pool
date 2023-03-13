
nsim    = 100000;
nburn   = round(0.4*nsim);

 horizon=2;
 
 path('..\DSGE_SW\results\sw2011',path);     
    load( strcat('sw2011pred_den_', num2str(horizon), '.mat'))
       PredDen_SW = pred_den(:,1);
       
[Tobs,n]    = size(PredDen_SW); 
T0_Forecast = horizon+1;
h_Forecast  = 0;

%% horizon = 2
 load('./results/MS_sample_02-Sep-2018_2.mat') 
 for k = 1:1:Tobs 
%     sort_filter(:,k) = sort(save_St(nburn+1:end,k),1);  
    sort_weight(:,k) = sort(save_Wt(nburn+1:end,k),1);  
 end
YY_mean_2     = mean(save_St(nsim-nburn:end,:),1)';
WW_mean_2     = mean(save_Wt(nsim-nburn:end,:),1)'; 
WW_median_2   = sort_weight(round((nsim-nburn)*0.5),:)'; 

%% horizon = 4
 load('./results/MS_sample_02-Sep-2018_4.mat') 
 for k = 1:1:Tobs 
%     sort_filter(:,k) = sort(save_St(nburn+1:end,k),1);  
    sort_weight(:,k) = sort(save_Wt(nburn+1:end,k),1);  
end
WW_mean_4     = mean(save_Wt(nsim-nburn:end,:),1)'; 
WW_median_4   = sort_weight(round((nsim-nburn)*0.5),:)'; 
YY_mean_4     = mean(save_St(nsim-nburn:end,:),1)';

%% horizon = 8
 load('./results/MS_sample_02-Sep-2018_8.mat') 
 for k = 1:1:Tobs 
%     sort_filter(:,k) = sort(save_St(nburn+1:end,k),1);  
    sort_weight(:,k) = sort(save_Wt(nburn+1:end,k),1);  
end
WW_mean_8     = mean(save_Wt(nsim-nburn:end,:),1)'; 
WW_median_8   = sort_weight(round((nsim-nburn)*0.5),:)'; 
YY_mean_8     = mean(save_St(nsim-nburn:end,:),1)';

%%

 varforecast = 6;
k = varforecast+1; 
disp(' Start Markov Switching Prediction Pool')

 ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs-h_Forecast)/4; 
%%
 
BC = csvread('./data/BC_Japan.csv',1,1);  

h_MS3 = figure('Position',[20,20,900,600],'File','Fig_5_MS',...
               'Name','Prob of Regime 1','Color','w');
title_name={ '(a) Probabilities of Regime that FA model is better','FontSize',14};
subplot(2,1,1);
hold on 
% hh = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[1.0,0.75,1.0]) ;
hh = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[0.5,0.65,1.0]) ;
      
ylim([0.0 1.0]);  
%    l1=plot(ti, YY_median(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','--','Color','r', 'LineWidth',2.5);
l1=plot(ti, YY_mean_2(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','--','Color','r', 'LineWidth',2.5);
l2=plot(ti, [nan(2,1);YY_mean_4(T0_Forecast+2:Tobs-h_Forecast,1)],'LineStyle','-','Color','k', 'LineWidth',2.5);
l3=plot(ti, [nan(6,1);YY_mean_8(T0_Forecast+6:Tobs-h_Forecast,1)],'LineStyle',':','Color','b', 'LineWidth',2.5);
title( title_name(1),'FontSize',14 )  
 legend([l1,l2,l3],{'2Q-Ahead','4Q-Ahead','8Q-Ahead'},'FontSize',11 );
hold off

subplot(2,1,2);
hold on 
% hh_w = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[1.0,0.75,1.0]) ;
hh_w = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[0.5,0.65,1]) ;
         
 l1=plot(ti, WW_mean_2(T0_Forecast:Tobs-h_Forecast,1), 'LineStyle','--','Color','r', 'LineWidth',2.5);   
 l2=plot(ti,[nan(2,1); WW_mean_4(T0_Forecast+2:Tobs-h_Forecast,1)], 'LineStyle','-','Color','k', 'LineWidth',2.5); 
 l3=plot(ti, [nan(6,1);WW_mean_8(T0_Forecast+6:Tobs-h_Forecast,1)], 'LineStyle',':','Color','b', 'LineWidth',2.5);
ylim([0.0 1.0]); 
title( '(b) Posterior Model Weight on FA Model','FontSize',14 );
legend([l1,l2,l3],{'2Q-Ahead','4Q-Ahead','8Q-Ahead'},'FontSize',11 );
hold off
 
%%
 est_date = datestr(date);   
name = ['./results/MS_Pool_state_',num2str(nsim),'_',est_date];
        saveas(h_MS3,name,'fig')       
          
save( ['./results/MS_sample_', est_date ,'.mat'], 'lambdasim', 'probsim',...
                                              'YY_median', 'save_St', 'save_Wt'  );  

%%