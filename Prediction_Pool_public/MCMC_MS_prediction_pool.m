% ------------------------------------------------------------------------
%   Markov Switching Predition Pooling 
% ------------------------------------------------------------------------
 clear all;
 clc;

 horizon=4;
 
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
pd          = 'Normal';
x0_old      = 0.25;                          %  x0_old = 0.2;
x1_old      = 0.75;                          %  x1_old = 0.5;
lambda0_old = x0_old ;
lambda1_old = x1_old ;
prob_old    = [0.95; 0.95];                  %  lambda_old = cdf(pd,x0_old,0,1);
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

save_St     = zeros(nsim,Tobs);
save_Wt     = zeros(nsim,Tobs);

 pmean =[0.05, 0.95]; % lambda0 lambda1
 pstdd  =[0.10, 0.10]; % standard deviation            
 pshape =[2,2]; % Gamma=2, Normal=3, Inverse Gamma=4 
% 
for j = 1:nsim
    
    % Step 1 - MH algorithm
      chk = 0;
    
    while chk == 0
      x0_new = x0_old + 0.05*randn(1,1);  
      x1_new = x1_old + 0.05*randn(1,1);      
      lambda0_new = x0_new;
      lambda1_new = x1_new;  
      
      if (x0_new < x1_new)
        chk = 1;
      end  
    end
    
    %{
    %x0_new + x1_new;  
    %lambda0_new = cdf(pd,x0_new,0,1);     
    %lambda1_new = cdf(pd,x1_new,0,1);
    %}
    
    sum_den =0;
    if (x0_new > 0) && (x0_new < 0.8)
        if  (x1_new > 0.2) && (x1_new < 1.0)
            for i = T0_Forecast:Tobs-h_Forecast
                if St_old(i)==0
                    den = (1 - lambda0_new)*PredDen_SW(i) + lambda0_new*PredDen_KK(i);
                else
                    den = (1 - lambda1_new)*PredDen_SW(i) + lambda1_new*PredDen_KK(i);
                end
                sum_den = sum_den + log(den);
            end            
            likenew = sum_den; 
        else
            likenew = -10^(1000);
        end
    else
        likenew = -10^(1000);
    end  
    
    para_new = [ lambda0_new; lambda1_new ];
     priornew =priodens(para_new, pmean, pstdd, pshape);
   
    para_old = [ lambda0_old; lambda1_old ];
     priorold =priodens(para_old, pmean, pstdd, pshape);
    
    % acceptance rate 
    postnew = likenew + priornew;
    postold = likeold + priorold;
    r = min(1,exp(postnew-postold));   
    if (rand < r)     
        likeold     =  likenew; 
        postold     =  postnew;             
        x0_old      =  x0_new;   
        x1_old      =  x1_new;
        lambda0_old =  lambda0_new;
        lambda1_old =  lambda1_new;
        acceptrate(j) = 1;            
    else
        acceptrate(j) = 0;
    end            
    lambdasim(j,:) = [lambda0_old lambda1_old] ;
         
    % Step 2 - Single move of Regimes (States)  
%     sampling_States_by_single_move;
   [St_old] = HamiltonFilter( PredDen_SW, PredDen_KK, lambda0_old, lambda1_old, ...
              prob_old,T0_Forecast,h_Forecast, horizon);
    
    save_St(j,:) = St_old;
    save_Wt(j,:) = lambda0_old*(1-St_old)+lambda1_old*St_old;
     
    % Step 3 - Sampling Transtion Probabilities
    pr00 = 9; pr01 = 1; pr11 = 9; pr10 = 1;
    N_00 = 0; N_01 = 0; N_10 = 0; N_11 = 0;

    for i = T0_Forecast:Tobs-1-h_Forecast   
        if (St_old(i)==0)
            if (St_old(i+1) ==0)   
                N_00 = N_00 + 1;
            else
                N_01 = N_01 + 1;        
            end
        else
            if (St_old(i+1) ==0)  
                N_10 = N_10 + 1;
            else
                N_11 = N_11 + 1;        
            end
        end
    end
    prob_old(1) = betarnd(N_00+pr00, N_01+pr01);
    prob_old(2) = betarnd(N_11+pr11, N_10+pr10); 
     probsim(j,:) = [ prob_old ] ;
     
   %==     display iteration number ===============     
     itr = num2str(j) ;      
     if mod(j, 2000)==0 
        disp([ itr 'th  iteration.......']);
     end
   %===============================================       
     
end

%% ------------------------------------------------------------------------
% Drawing figures
% ------------------------------------------------------------------------ 
   
h_MS1 = figure('Position',[20,20,900,600],'Name','Posterior Density of Lambda_1 and Lambda_2','Color','w');
subplot(2,2,1)
[density,x1]  = ksdensity(lambdasim(nburn:end,1));
plot(x1,density,'LineStyle','-','Color','r','LineWidth',2.5);
title( 'Posterior Density of \lambda_1 under Regime 1','FontSize',12 );
xlim([0 1]);

subplot(2,2,2)
[density,x1]  = ksdensity( lambdasim(nburn:end,2) );
plot(x1,density,'LineStyle','-','Color','r','LineWidth',2.5);
title( 'Posterior Density of \lambda_2 under Regime 2','FontSize',12 );
xlim([0 1]);

subplot(2,2,3)
plot(lambdasim(nburn:end,1));
title( 'Trace of \lambda_1 under Regime 1','FontSize',12 );
ylim([0 1]);

subplot(2,2,4)
plot( lambdasim(nburn:end,2) );
title( 'Trace of \lambda_2  under Regime 2','FontSize',12 );
ylim([0 1]);

%{
% figure('Name','Prob of two Regimes ','Color','w');
% subplot(2,1,1)
% plot(ti, mean(save_St(:,T0:Tobs-h),1))
% title( 'Prob of Regime 0','FontSize',14 );
% subplot(2,1,2)
% plot(ti, 1-mean(save_St(:,T0:Tobs-h),1))
% title( 'Prob of Regime 1','FontSize',14 );
%}

lambda0 = mean(lambdasim(nburn:end,1));
lambda1 = mean(lambdasim(nburn:end,2)); % + mean(lambdasim(nburn:end,1));
     
%% *********************************************
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

%% ------------------------------------------------------------------------
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

%% ------------------------------------------------------------------------  
%  drawing graph  
% ------------------------------------------------------------------------  

ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs-h_Forecast)/4; 

h_MS2 = figure('Position',[20,20,900,600],'Name','Prob of Regime 1','Color','w');
title_name={ '(a) Probabilities of Regime that FF model is better, \lambda_2','FontSize',14};
subplot(2,1,1);
hold on 
hh = area(ti,  YY_f_save(T0_Forecast:Tobs-h_Forecast,:)) ; 
      
for j = 1:m    
    set(hh(j),'FaceColor',...
    [col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3) ])    
    set(hh(2*(m+1)-j),'FaceColor',...
    [col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3) ])
end
set(hh(m+1),'FaceColor',[1-col1 1-col2 1-col3 ])
set(hh(1),'FaceColor',[1 1 1 ])    
set(hh,'LineStyle','none') % Set all to same value
      
ylim([0.0 1.0]);  
   plot(ti, YY_median(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','--','Color','r', 'LineWidth',2.5);
   plot(ti, YY_mean(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','k', 'LineWidth',2.5);
title( title_name(1),'FontSize',14 )  
hold off

subplot(2,1,2);
hold on 
hh_w = area(ti, WW_save(T0_Forecast:Tobs-h_Forecast,:)) ; 
      
for j = 1:m    
    set(hh_w(j),'FaceColor',...
    [col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3) ])    
    set(hh_w(2*(m+1)-j),'FaceColor',...
    [col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3) ])
end
set(hh_w(m+1),'FaceColor',[1-col1 1-col2 1-col3 ])
set(hh_w(1),'FaceColor',[1 1 1 ])    
set(hh_w,'LineStyle','none')

%lambda_mean =  lambda0*(1-YY_mean(T0:Tobs-h,1))+lambda1*(YY_mean(T0:Tobs-h,1));
          
 plot(ti, WW_median(T0_Forecast:Tobs-h_Forecast,1), 'LineStyle','--','Color','r', 'LineWidth',2.5);   
 plot(ti, WW_mean(T0_Forecast:Tobs-h_Forecast,1), 'LineStyle','-','Color','k', 'LineWidth',2.5); 
ylim([0.0 1.0]); 
title( '(b) Posterior Model Weight on FF Model, \lambda_t','FontSize',14 );
hold off

 save( './data/MS_WW_mean.mat', 'WW_mean','WW_median')  ; 
 
 
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
   l1=plot(ti, YY_median(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','--','Color','r', 'LineWidth',2.5);
   l2=plot(ti, YY_mean(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','k', 'LineWidth',2.5);
title( title_name(1),'FontSize',14 )  
legend([l1,l2],{'Median','Mean'},'FontSize',11 );
hold off

subplot(2,1,2);
hold on 
% hh_w = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[1.0,0.75,1.0]) ;
hh_w = area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(100),'LineStyle','non','FaceColor',[0.5,0.65,1]) ;
         
 l1=plot(ti, WW_median(T0_Forecast:Tobs-h_Forecast,1), 'LineStyle','--','Color','r', 'LineWidth',2.5);   
 l2=plot(ti, WW_mean(T0_Forecast:Tobs-h_Forecast,1), 'LineStyle','-','Color','k', 'LineWidth',2.5); 
ylim([0.0 1.0]); 
title( '(b) Posterior Model Weight on FA Model','FontSize',14 );
legend([l1,l2],{'Median','Mean'},'FontSize',11 );
hold off
 
%%
 est_date = datestr(date);   
name = ['./results/MS_Pool_state_',num2str(nsim),'_',est_date,'_',num2str(horizon)];
        saveas(h_MS3,name,'fig')       
          
save( ['./results/MS_sample_', est_date ,'_',num2str(horizon),'.mat'], 'lambdasim', 'probsim',...
                                              'YY_median', 'save_St', 'save_Wt'  );  

%%
%----------------------------------------------------------------
%   Summary of Posterior Estimates of Lambda
%-----------------------------------------------------------------     
     mean_lam1 = mean(lambdasim(nburn+1:end,1));
     std_lam1  = std(lambdasim(nburn+1:end,1));
     
     mean_lam2 = mean(lambdasim(nburn+1:end,2));
     std_lam2  = std(lambdasim(nburn+1:end,2));
     
     mean_prob1 = mean(probsim(nburn+1:end,1));
     std_prob1  = std(probsim(nburn+1:end,1));
     
     mean_prob2 = mean(probsim(nburn+1:end,2));
     std_prob2  = std(probsim(nburn+1:end,2));
     
     para_save = zeros(nsim-nburn, 4 );
     
     for i = 1:1:4
        sort_para=zeros( nsim-nburn, 1 );
        if i < 3
          sort_para(:,1) = sort(lambdasim(nburn+1:end,i),1);
        else
          sort_para(:,1) = sort(probsim(nburn+1:end,i-2),1);
        end  
        
          a  =0.80;   % (1-2*rate) percentage of Credibile interval of Parameters  
          rate = (1-a)/2;
          para_low(i) = sort_para(round((nsim-nburn)*rate),1); 
          para_up(i)  = sort_para(round((nsim-nburn)*(1-rate)),1);
    
        if i < 3
           para_save(:,i) = lambdasim(nburn+1:end,i);
        else
           para_save(:,i) = probsim(nburn+1:end,i-2);
        end
        
          p=1; n=1;
          
           cal_inefficiency_MS    
%           inefficiency(i)=0;
    
     end
     
     disp('  mean  s.d. [ lower Band Upper Band ] inefficiency');
     disp([ mean_lam1 std_lam1 para_low(1) para_up(1) inefficiency(1) ] );
     disp([ mean_lam2 std_lam2 para_low(2) para_up(2) inefficiency(2) ] );
     disp([ mean_prob1 std_prob1 para_low(3) para_up(3) inefficiency(3) ] );
     disp([ mean_prob2 std_prob2 para_low(4) para_up(4) inefficiency(4) ] );
     
%% -------------------------------------------------------------
%    output to text file
%-------------------------------------------------------------

para_name_MS = char('\lambda0', '\lambda1', 'p11', 'p22'); 

est_date = datestr(date);   
result_name = ['./results/MS_pool', est_date ,'_',num2str(horizon), '.txt'];          
fileID = fopen(result_name,'w');
   
fprintf(fileID,'\n---------------------------------------------------------------');
fprintf(fileID,'\n\n                        [ESTIMATION RESULT]');
fprintf(fileID,'\n---------------------------------------------------------------');

fprintf(fileID, '\n Parameter   mean      s.d.     [lower     Upper] inefficiency \n');
fprintf(fileID, '%s %9.3f %9.3f %9.3f %9.3f %9.3f  \n', ...
        para_name_MS(1,:), ... 
        [ mean_lam1 std_lam1 para_low(1) para_up(1) inefficiency(1) ] );
fprintf(fileID, '%s %9.3f %9.3f %9.3f %9.3f %9.3f  \n', ...
        para_name_MS(2,:), ... 
        [ mean_lam2 std_lam2 para_low(2) para_up(2) inefficiency(2) ]  );   
fprintf(fileID, '%s %9.3f %9.3f %9.3f %9.3f %9.3f  \n', ...
        para_name_MS(3,:), ... 
       [ mean_prob1 std_prob1 para_low(3) para_up(3) inefficiency(3) ] );
fprintf(fileID, '%s %9.3f %9.3f %9.3f %9.3f %9.3f  \n', ...
       para_name_MS(4,:), ... 
      [ mean_prob2 std_prob2 para_low(4) para_up(4) inefficiency(4) ]  );    

fprintf(fileID,'\n --------------------------------------------------------------');
fprintf(fileID,'\n --------------------------------------------------------------');
fclose(fileID);

est_date = datestr(date);   
name = ['./results/MS_Pool_para',num2str(nsim),'_',est_date,'_',num2str(horizon)];
        saveas(h_MS1,name,'fig')
        


