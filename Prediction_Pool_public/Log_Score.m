 clear;
%  clc;

horizon=2;

varforecast = 6;
k = varforecast+1;

%   period = 'bubble'
period = 'full_sample'
% period = 'Pre_Bubble'
% period = 'Post_Bubble'

 % output=1, % cons = 2   % inv=3   % real wage=4   %inflation=5
title_var = { 'output', 'cons', 'inv', 'real wage', 'inflation'};

disp(' Start Calculating Log Score of Prediction Pool Methods')

path('function',path);
path('..\DSGE_KK2\results\FF2011',path);
    load(strcat('FF2011pred_den_', num2str(horizon), '.mat'))
       PredDen_KK = pred_den(:,k);
path('..\DSGE_SW\results\sw2011',path);     
    load( strcat('sw2011pred_den_', num2str(horizon), '.mat'))
       PredDen_SW = pred_den(:,k);
       
 disp('Loading new prediction density');
    
   
%     load 'OP2011pred_den.mat'
%        PredDen_OP = pred_den;

  BC = csvread('./data/BC_Japan.csv',1,1);   

[Tobs,n] = size(PredDen_SW); 

switch period
    case 'full_sample'        
         T0_Forecast = horizon+1;
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

% ------------------------------------------------------------------------
%   i.) Log Score of Static Pooling
% ------------------------------------------------------------------------
if horizon == 8
  load('./data/mean_lam_8.mat')
elseif horizon == 4
  load('./data/mean_lam_4.mat')
elseif horizon == 2
  load('./data/mean_lam_2.mat')
end  
lambda_Static   = mean_lam;  
switch period
    case 'full_sample'        
        lambda_Static   = mean_lam;
    case 'bubble'
         lambda_Static   = 0.576;
    case 'Pre_Bubble'
         lambda_Static   = 0.773;
    case 'Post_Bubble'
         lambda_Static   =  0.313;
end         
Pred_Static     = zeros(Tobs,1);

log_score_Static    = 0;
log_score_SW      = 0;
log_score_KK      = 0;

for i = T0_Forecast:Tobs-h_Forecast
    Pred_Static(i)      =   (1 - lambda_Static)*PredDen_SW(i) + lambda_Static*PredDen_KK(i);
    log_score_SW      =   log_score_SW + log(PredDen_SW(i));
    log_score_KK      =   log_score_KK + log(PredDen_KK(i));
    log_score_Static    =   log_score_Static + log( Pred_Static(i));  
end   
disp( ['log_score_SW      = ' num2str(log_score_SW)] );
disp( ['log_score_KK      = ' num2str(log_score_KK)] );
disp( ['log_score_Static    = ' num2str(log_score_Static)] );

% ------------------------------------------------------------------------
%   ii.) Log Score of MS Pooling
% ------------------------------------------------------------------------
if horizon == 8
  load('./results/MS_sample_02-Sep-2018_8.mat')
elseif horizon == 4
  load('./results/MS_sample_02-Sep-2018_4.mat')
elseif horizon == 2
  load('./results/MS_sample_02-Sep-2018_2.mat')
end  
  
Pred_MS         =   zeros(Tobs,1);
lambda_MS       =   mean(save_Wt,1);
%  lambda_MS       =   WW_mean;
log_score_MS    =   0;
for i = T0_Forecast:Tobs-h_Forecast
    Pred_MS(i)      = (1 - lambda_MS(i))*PredDen_SW(i) + lambda_MS(i)*PredDen_KK(i);
    log_score_MS    = log_score_MS + log(Pred_MS(i));
end   
disp( ['log_score_MS = ' num2str(log_score_MS)] );


% % ------------------------------------------------------------------------
% %  iii.) Log Score of Dynamic Pooling
% % ------------------------------------------------------------------------
% load './data/DP_WW_mean2.mat'
% 
% Pred_DP2         =   zeros(Tobs,1);
% lambda_DP2       =   YY_median;
% %  lambda_DP2       =   YY_mean;
% log_score_DP2    =   0;
% 
% for i = T0_Forecast:Tobs-h_Forecast
%     Pred_DP2(i) = (1-lambda_DP2(i))*PredDen_SW(i)+lambda_DP2(i)*PredDen_KK(i);
%     log_score_DP2 = log_score_DP2 + log( Pred_DP2(i));
% end   
% disp( ['log_score_DP_1_para = ' num2str(log_score_DP2)] );
% 

% ------------------------------------------------------------------------
%  iii.) Log Score of Dynamic Pooling
% ------------------------------------------------------------------------
if horizon == 8
  load('./data/DP_WW_mean1_8.mat')
elseif horizon == 4
 load('./data/DP_WW_mean1_4.mat')
elseif horizon == 2
  load('./data/DP_WW_mean1_2.mat')
end  

Pred_DP         =   zeros(Tobs,1);
 lambda_DP       =   YY_mean;
%  lambda_DP       =   YY_median;
log_score_DP    =   0;

for i = T0_Forecast:Tobs-h_Forecast
    Pred_DP(i) = (1-lambda_DP(i))*PredDen_SW(i)+lambda_DP(i)*PredDen_KK(i);
    log_score_DP = log_score_DP + log( Pred_DP(i));
end   
disp( ['log_score_DP_3_para = ' num2str(log_score_DP)] );


% 
% 
 ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs-h_Forecast)/4;
% 
% h_log_score1 = figure('Position',[20,20,900,600],'Name','Log Score of SW, KK and 3 methods'); hold on
%   area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(-100),'LineStyle','non','FaceColor',[1.0,0.75,1.0]) 
% plot(ti, log(PredDen_SW(T0_Forecast:Tobs-h_Forecast)), 'LineStyle',':','Color','r','LineWidth',1.5);
% plot(ti, log(PredDen_KK(T0_Forecast:Tobs-h_Forecast)), 'LineStyle',':','Color','b','LineWidth',1.5);      
% plot(ti, log(Pred_Static(T0_Forecast:Tobs-h_Forecast)), 'LineStyle','--','Color','b','LineWidth',1.5);
% plot(ti, log(Pred_MS(T0_Forecast:Tobs-h_Forecast)), 'LineStyle','-','Color','k','LineWidth',1.5);
% plot(ti, log(Pred_DP(T0_Forecast:Tobs-h_Forecast)), 'LineStyle','-','Color','g','LineWidth',1.5);
%  ylim([-18,-4])
% title( 'Log Score of NK, FF, and 3 pooling methods','FontSize',12 );
% legend('Recessions','NK', 'FF', 'Static method', 'MS method','Dynamic method' );
% 
% hold off    
% 
h_log_score2 =figure('Position',[20,20,900,300],'Name','Predictive Densities of NK model and FF model'); hold on
%  area(ti', BC(T0_Forecast:end)*(-100),'LineStyle','non','FaceColor',[0,0.90,0.90]) 
  area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(-100),'LineStyle','non','FaceColor',[1.0,0.75,1.0]) 
 plot(ti, log(PredDen_KK(T0_Forecast:Tobs-h_Forecast)), 'LineStyle','-','Color','b','LineWidth',1.5);  
plot(ti, log(PredDen_SW(T0_Forecast:Tobs-h_Forecast)), 'LineStyle','--','Color','r','LineWidth',1.5);     
% title( [ 'Predictive Densities of NK model with and without Financial Friction'],'FontSize',12 );
legend({'Recessions','FA model','NK model'},'FontSize',12  );
  ylim([-18,-4])
  xlim([1980, 2000])
hold off    
% 
% 
% %    plot(ti, log(Pred_Static(T0:Tobs-h)), 'LineStyle','--','Color','b','LineWidth',1.5);
% %    plot(ti, log(Pred_MS(T0:Tobs-h)), 'LineStyle','--','Color','r','LineWidth',1.5);
% %    plot(ti, log(Pred_DP(T0:Tobs-h)), 'LineStyle','-','Color','black','LineWidth',1.5);
% %    plot(ti, [ log(PredDen_SW(T0:Tobs-h)) log(PredDen_KK(T0:Tobs-h))...
% %              log(Pred_MS(T0:Tobs-h)) log(Pred_DP(T0:Tobs-h))] );


%-------------------------------------------------------------
%    output to text file
%-------------------------------------------------------------

est_date = datestr(date);   
result_name = ['c:\DSGE\DSGE_KK/results/Log_score_', est_date , '.txt'];          
fileID = fopen(result_name,'w');
   
fprintf(fileID,'\n---------------------------------------------------------------');
fprintf(fileID,'\n\n                        [ESTIMATION RESULT]');
fprintf(fileID,'\n---------------------------------------------------------------');

fprintf(fileID, ['\n log_score_NK      = ' num2str(log_score_SW) ] );
fprintf(fileID, ['\n log_score_FF      = ' num2str(log_score_KK)] );
fprintf(fileID, ['\n log_score_Static  = ' num2str(log_score_Static)] );
fprintf(fileID, ['\n log_score_MS      = ' num2str(log_score_MS)] );
fprintf(fileID, ['\n log_score_DP1      = ' num2str(log_score_DP)] );
% fprintf(fileID, ['\n log_score_DP2      = ' num2str(log_score_DP2)] );

fprintf(fileID,'\n -----------------------------------');
fprintf(fileID,'\n -----------------------------------');
fclose(fileID);

est_date = datestr(date);   
% name = ['./results/Log_score_1_',est_date];
%         saveas(h_log_score1,name,'fig')
%         
% est_date = datestr(date);   
% name = ['./results/Log_score_2_',est_date];
%         saveas(h_log_score2,name,'fig')   


%% End
