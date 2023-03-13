


clear;
%  clc;

horizon=8;

k = 1; 
varforecast = 6;
 % output=1, % cons = 2   % inv=3   % real wage=4   %inflation=5
title_var = { 'output', 'cons', 'inv', 'real wage', 'inflation'};

disp(' Start Calculating Log Score of Prediction Pool Methods')

path('function',path);
path('..\DSGE_KK2\results\FF2011',path);
    load(strcat('FF2011pred_den_', num2str(horizon), '.mat'))
       PredDen_KK = pred_den;
 path('..\DSGE_SW\results\sw2011',path);     
    load( strcat('sw2011pred_den_', num2str(horizon), '.mat'))
       PredDen_SW = pred_den;
       
 display('Loading new prediction density');  


  BC = csvread('./data/BC_Japan.csv',1,1);   

[Tobs,n] = size(PredDen_SW); 
T0_Forecast = 5;
h_Forecast  = 0;

ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs)/4;

h_log_score1 =figure('Position',[20,20,900,600],'Name','Predictive Densities of NK model and FF model');

for i =1:varforecast
subplot(varforecast,1,i)
 hold on
%  area(ti', BC(T0_Forecast:end)*(-100),'LineStyle','non','FaceColor',[0,0.90,0.90]) 
 area(ti', BC(T0_Forecast:end)*(-100),'LineStyle','non','FaceColor',[1.0,0.75,1.0]) 
plot(ti, log(PredDen_KK(T0_Forecast:Tobs-h_Forecast,i)), 'LineStyle','-','Color','b','LineWidth',1.5);      
plot(ti, log(PredDen_SW(T0_Forecast:Tobs-h_Forecast,i)), 'LineStyle','--','Color','r','LineWidth',1.5);


if i == 1   ylim([-5,-0]); title( [ 'Output'],'FontSize',12 );
elseif i ==2   ylim([-5,-0]);  title( [ 'Consumption'],'FontSize',12 );
elseif i ==3   ylim([-10,-0]); title( [ 'Investment'],'FontSize',12 );
elseif i ==4    ylim([-4,-0]); title( [ 'Real Wage'],'FontSize',12 );
elseif i ==5      ylim([-1.2,-0]); title( [ 'Inflation'],'FontSize',12 );
elseif i ==6      ylim([-1,-0]); title( [ 'Nominal Rate'],'FontSize',12 );    
    legend('Recessions', 'FA model','NK model'  );
end    
hold off  
end

est_date = datestr(date);   
name = ['c:\DSGE\DSGE_KK/results/Log_score_sngl_',est_date];
        saveas(h_log_score1,name,'fig')
       
        
        load 'FF2011pred_diff.mat'
       Pred_diff_KK = save_diff;
load 'sw2011pred_diff.mat'
       Pred_diff_SW = save_diff;
       
var_title ={'output','cons','inv','wage','inf','Nominal rate' };
h_log_score3 =figure('Position',[20,20,900,600],'Name','Predictive Densities of NK model and FF model');
for i = 1:6
  subplot(6,1,i)
  hold on
    hh=area(ti',[-100*ones(size(ti,2),1) BC(T0_Forecast:end)*(200)]);
     set(hh(1),'FaceColor',[1 1 1 ]);
     set(hh(2),'FaceColor',[1 0.75 1 ]);
     set(hh,'LineStyle','none');
     
     plot(ti,Pred_diff_KK(T0_Forecast:Tobs-h_Forecast,i),'LineStyle','-','Color','b','LineWidth',1.5 )
    plot(ti,Pred_diff_SW(T0_Forecast:Tobs-h_Forecast,i),'LineStyle','--','Color','r','LineWidth',1.5 )
   
%     plot(ti,zeros(size(ti,2),1),'r--' )
  hold off 
  if i == 1   ylim([-3.0,2.5]); title( [ 'Output'],'FontSize',12 );
elseif i ==2   ylim([-3.0,2.5]);  title( [ 'Consumption'],'FontSize',12 );
elseif i ==3   ylim([-7.5,7.5]); title( [ 'Investment'],'FontSize',12 );
elseif i ==4    ylim([-2.0,1.5]); title( [ 'Real Wage'],'FontSize',12 );
elseif i ==5      ylim([-0.8,0.8]); title( [ 'Inflation'],'FontSize',12 );
elseif i ==6      ylim([-0.5,0.5]); title( [ 'Nominal Rate'],'FontSize',12 );    
    legend('','Recessions','FA model','NK model'  );
end    
end   

est_date = datestr(date);   
name = ['c:\DSGE\DSGE_KK/results/Fore_Err_sngl_',est_date];
        saveas(h_log_score3,name,'fig')

% var_title ={'output','cons','inv','wage','inf','Nominal rate' };
% h_log_score4 =figure('Position',[20,20,900,600],'Name','Predictive Densities of NK model and FF model');
% for i = 4:6
%   subplot(3,1,i-3)
%   hold on
%     area(ti', BC(T0_Forecast:end)*(-100),'LineStyle','non','FaceColor',[0,0.90,0.90]) 
%      plot(ti,Pred_diff_KK(T0_Forecast:Tobs-h_Forecast,i),'LineStyle','-','Color','b','LineWidth',1.5 )
%     plot(ti,Pred_diff_SW(T0_Forecast:Tobs-h_Forecast,i),'LineStyle','--','Color','r','LineWidth',1.5 )
%    
% %     plot(ti,zeros(size(ti,2),1),'r--' )
%   hold off 
%   if i == 1   ylim([-3,3]); title( [ 'Output'],'FontSize',12 );
% elseif i ==2   ylim([-5,5]);  title( [ 'Consumption'],'FontSize',12 );
% elseif i ==3   ylim([-10,10]); title( [ 'Investment'],'FontSize',12 );
% elseif i ==4    ylim([-2,2]); title( [ 'Real Wage'],'FontSize',12 );
% elseif i ==5      ylim([-2,2]); title( [ 'Inflation'],'FontSize',12 );
% elseif i ==6      ylim([-2,1]); title( [ 'Nominal Rate'],'FontSize',12 );    
%     legend('Recessions','NK model', 'Financial Friction','FontSize',12  );
% end    
% end    

%  period = 'bubble'
% period = 'full_sample'
% period = 'Pre_Bubble'
period = 'Post_Bubble'

switch period
    case 'full_sample'        
         T0_Forecast = 5;
         h_Forecast  = 0;
    case 'bubble'
         T0_Forecast = 25;         
         h_Forecast  = Tobs - 51;
    case 'Pre_Bubble'
         T0_Forecast = 5;         
         h_Forecast  = Tobs - 23;   
    case 'Post_Bubble'
         T0_Forecast = 52;         
         h_Forecast  = 0; 
end

num_b = 0; 
num_r = 0; 
boom1 = zeros(6,1);
recession1 = zeros(6,1);
boom2 = zeros(6,1);
recession2 = zeros(6,1);

for i = T0_Forecast:Tobs-h_Forecast
  if BC(i)==0
       boom1 = boom1 + Pred_diff_KK(i,:)';
       boom2 = boom2 + Pred_diff_SW(i,:)';
       num_b = num_b + 1;
  else
       recession1 = recession1 + Pred_diff_KK(i,:)';
       recession2 = recession2 + Pred_diff_SW(i,:)';
       num_r = num_r + 1;
  end
end  
  
disp('FF model')
disp( [ (boom1/num_b)' ] )
disp( [ (recession1/num_r)' ] )

disp('NK model')
disp( [ (boom2/num_b)' ] )
disp( [ (recession2/num_r)' ] )


