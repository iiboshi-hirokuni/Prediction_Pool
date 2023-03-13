%% ------------------------------------------------------------------------
%%  Dynamic Prediction Pooling 
%% ------------------------------------------------------------------------
%    clear all;
%    clc;
% % 
% 
% horizon=8;
% 
%  tic
%  varforecast = 6;
% k = varforecast+1; 
% %% 並列処理の　オン・オフ
% %  parallel = 'off' 
% parallel = 'on' 
% ncores = 4;
% %        delete(gcp('nocreate')) %   
% %         parpool('local', ncores)
% %%
% %   parameter = 'all'    %% estimate three parameters: rho, mu, sigma 
% parameter = 'one'   %% estimate one parameters: only rho
% 
% disp(' Start Dynamic Prediction Pool')
%  disp('Loading new prediction density');
% path('function',path);
% path('..\DSGE_KK2\results\FF2011',path);
%     load(strcat('FF2011pred_den_', num2str(horizon), '.mat'))
%        PredDen_KK = pred_den(:,k);
%   path('..\DSGE_SW\results\sw2011',path);     
%     load( strcat('sw2011pred_den_', num2str(horizon), '.mat'))
%        PredDen_SW = pred_den(:,k);
% %     load 'OP2011pred_den.mat'
% %        PredDen_OP = pred_den;
% 
% [Tobs,n]   = size(PredDen_KK); 
% T0         = horizon;
% nparticle   = 500;
%  nsim        = 500;
%  start_sim   = 1;
%  nburn       = ceil(0.3*nsim);
% rho_old     = 0.5;
% mu_old      = 0; % icdf('norm',0.6,0,1);
% sigma_old   =  1.0;
% T0_Forecast = horizon+1;
% h_Forecast  = 0;
% 
% switch parameter
% case  'all' 
%  rho_old     = 0.5;  
%  sigma_old   =  1.5;   
%   c1 = 0.05;
%   c2 = 0.01;
%   c3 = 0.05;
% case 'one'
%  rho_old     = 0.75;  
%   c1 = 0.05;
%   c2 = 0;
%   c3 = 0;
% end
% 
% %% prior setting
% switch parameter
% case  'all' 
%  pmean =[0.5, 0.5, 1 ]; % rho, mu, sigma 
%  pstdd  =[0.15, 1, 10]; % standard deviation            
%  pshape =[1,3.4]; % beta=1, Normal=3, Inverse Gamma=4 
% case  'one' 
%  pmean =[0, 0, 0 ]; % rho, mu, sigma 
%  pstdd  =[1, 0, 0]; % standard deviation            
%  pshape =[5,0.0]; % beta=1, Normal=3, Inverse Gamma=4, uniform=5 
% end
% 
% likeold     = -10^(1000);
% postold     = -10^(1000);
% priornew    = 0;
% priorold    = 0;
% 
%  lambdasim = zeros(nsim,4);
%  save_filter = zeros( nsim, Tobs );
% 
% % switch parallel  
% %         case 'on'    
% %         poolobj = parpool;
% % end
% 
% 
% 
% rv = randn(Tobs,1,nparticle );

%%
   start_sim = 850;
  rho_old = lambdasim(start_sim,1);
  mu_old = lambdasim(start_sim,2);
  sigma_old = lambdasim(start_sim,3);
  
  %    nsim        = 5000;
%    nburn       = ceil(0.3*nsim);
% %   lambdasim = [lambdasim; zeros(nsim-start_sim,4)];
% %    save_filter =[save_filter; zeros( nsim-start_sim, Tobs )];
%%

for j = start_sim+1:nsim      
    rho_new  =  rho_old + c1*randn(1,1);
    mu_new   =  mu_old  + c2*randn(1,1);
    sigma_new = sigma_old + c3*randn(1,1);
        
%     sum_den     =   0;
    if (rho_new > 0)&&(rho_new < 1)&&(sigma_new>0)          
      
      switch parallel  
        case 'on'    
%         poolobj = parpool;
        %% Perticle Filter
         [resample_xt, log_lik, mean_xt] = ParticleFilter_par(PredDen_SW, PredDen_KK, T0_Forecast,...
                                           h_Forecast, nparticle, rho_new, mu_new, sigma_new, rv,ncores);
       case 'off'
         [resample_xt, log_lik, mean_xt] = ParticleFilter(PredDen_SW, PredDen_KK, T0_Forecast,...
                                           h_Forecast, nparticle, rho_new, mu_new, sigma_new, rv);
       end                             
      
      likenew = log_lik;                   
                            
%       lambda_DP       =   mean_xt;
%       log_score_DP    =   0;
%        for i = T0_Forecast:Tobs-h_Forecast
%            Pred_DP(i) = (1-lambda_DP(i))*PredDen_SW(i)+lambda_DP(i)*PredDen_KK(i);
%            log_score_DP = log_score_DP + log( Pred_DP(i));
%        end                       
%        likenew = log_score_DP;
     
    else
        likenew = -10^(1000);
    end 
    
    % --------------------------------------------------------------------
    %   MH algorithm - compute the acceptance rate 
    % --------------------------------------------------------------------
    para_new = [rho_new; mu_new; sigma_new];
     priornew =priodens(para_new, pmean, pstdd, pshape);
   
    para_old = [rho_old; mu_old; sigma_old];
     priorold =priodens(para_old, pmean, pstdd, pshape);
    
    postnew = likenew + priornew;      
    postold = likeold + priorold;       
    r = min(1,exp(postnew-postold));   
    if (rand < r)     
        likeold         =  likenew; 
        postold         =  postnew;             
        re_xt_old       = resample_xt(:,ceil(rand*nparticle));
        
        rho_old     =  rho_new;
        mu_old      =  mu_new;
        sigma_old   =  sigma_new;
        
        acceptrate(j)   =  1;       
    else
        acceptrate(j)   =   0;
    end
    
    lambdasim(j,:) = [ rho_old  cdf('norm',mu_old,0,1) sigma_old postold];
    save_filter(j,:)    = re_xt_old;
    
     %==     display iteration number ===============     
     itr = num2str(j) ;     
           
     if mod(j, nsim/100)==0 
        disp([ itr 'th  iteration.......']);
%         if j > nburn
          disp([ 'accept rate = ' num2str(mean(acceptrate(:)) ) ] );
          
          disp([ 'rho, Phi(mu), sigma = ' num2str(mean(lambdasim(ceil(j/2):j,:),1)) ] );
%         end
       
          disp( [ 'time = ' num2str(floor(toc/60/60)) ' (hr)' ...
                    num2str(floor(mod(toc/60, 60))) '(min)'...
                    num2str(floor(mod(toc,60))) '(sec)']);

     end
   %===============================================     
    
    
end    

% switch parallel  
%         case 'on'    
%          delete(poolobj);
%  end

%%

% nsim        = 950;
%  nburn       = ceil(0.3*nsim);
% lambdasim = lambdasim(1:5000,:);
% save_filter= save_filter(1:5000,:);
%%
% switch parameter
%     case 'all'  
%       save('./data/dynamic_1_sampling.mat', 'lambdasim', 'save_filter');
%     case 'one'
%       save('./data/dynamic_2_sampling.mat', 'lambdasim', 'save_filter');
% end      

%%

ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs)/4; 
   
figure('Position',[20,20,900,600],'Name','Static Optimal Pools (Weight on KK)','Color','w');
subplot(1,3,1)
[density,x1] = ksdensity(lambdasim(nburn:end,1));
plot(x1,density,'LineStyle','-','Color','b','LineWidth',2.5);
title( 'Posterior Density of \rho','FontSize',12 );
subplot(1,3,2)
[density,x1] = ksdensity(lambdasim(nburn:end,2));
plot(x1,density,'LineStyle','-','Color','b','LineWidth',2.5);
title( 'Posterior Density of \mu','FontSize',12 );
subplot(1,3,3)
[density,x1] = ksdensity(lambdasim(nburn:end,3));
plot(x1,density,'LineStyle','-','Color','b','LineWidth',2.5);
title( 'Posterior Density of \sigma','FontSize',12 );

% plot(lambdasim(nburn:end));
% title('Trace of \lambda','FontSize',12 );


% figure('Position',[20,20,900,600],'Name','Prediction Score of KK and SW','Color','w');
%  plot(ti', log(PredDen_SW(T0_Forecast:Tobs-h_Forecast)), 'r', 'LineWidth',2)
% hold on
%     plot(ti', log(PredDen_KK(T0_Forecast:Tobs-h_Forecast)), 'b--', 'LineWidth',2) 
% %     plot(ti', log(PredDen_OP(T0_Forecast:Tobs-h_Forecast)),'g-', 'LineWidth',2) 
% hold off
%         
% title( 'Predictive Densities of KK and SW','FontSize',12 );
% legend('SW', 'KK');

%--------------------------------------------------------------------------
%   Summary of Posterior Estimates of Lambda  
%--------------------------------------------------------------------------   

sort_para       = zeros(nsim-nburn, 3);
for i = 1:3
  mean_lam(i)        = mean(lambdasim(nburn+1:end,i));  
  std_lam(i)         = std(lambdasim(nburn+1:end,i),1);
  sort_para(:,i)  = sort(lambdasim(nburn+1:end,i),1);
end  
    
% Calculating of Posterior estimates 
a           =   0.90;   % (1-2*rate) percentage of Credibile interval of Parameters  
rate        =   (1-a)/2;
for i = 1:3
 para_low(i)    =   sort_para(ceil((nsim-nburn)*rate),i); 
 para_up(i)     =   sort_para(ceil((nsim-nburn)*(1-rate)),i);
end 

para_save   =   lambdasim(nburn:end,:);
p           =   3;
cal_inefficiency    
disp('  mean  s.d. [ lower Band Upper Band ] inefficiency');
disp([ mean_lam(1) std_lam(1) para_low(1) para_up(1) inefficiency(1) ] );
disp([ mean_lam(2) std_lam(2) para_low(2) para_up(2) inefficiency(2) ] );
disp([ mean_lam(3) std_lam(3) para_low(3) para_up(3) inefficiency(3) ] );

% save ./data/mean_lam.mat  mean_lam; 

%% -------------------------------------------------------------
%    output to text file
%-------------------------------------------------------------

para_name_Dy = char('\rho', '\mu', '\sigma'); 

est_date = datestr(date);   
switch parameter
    case 'all'    
     result_name = ['./results/Dynamic_pool1', est_date,'_', num2str(horizon) , '.txt'];          
    case 'one'    
     result_name = ['./results/Dynamic_pool2', est_date,'_', num2str(horizon) , '.txt']; 
end
    fileID = fopen(result_name,'w');
   
fprintf(fileID,'\n---------------------------------------------------------------');
fprintf(fileID,'\n\n                        [ESTIMATION RESULT]');
fprintf(fileID,'\n---------------------------------------------------------------');

fprintf(fileID, '\n Parameter   mean      s.d.     [lower     Upper] inefficiency \n');
fprintf(fileID, '%s %9.3f %9.3f %9.3f %9.3f %9.3f  \n', ...
        para_name_Dy(1,:), ... 
        [ mean_lam(1) std_lam(1) para_low(1) para_up(1) inefficiency(1) ] );
fprintf(fileID, '%s %9.3f %9.3f %9.3f %9.3f %9.3f  \n', ...
        para_name_Dy(2,:), ... 
        [ mean_lam(2) std_lam(2) para_low(2) para_up(2) inefficiency(2) ]  );   
 fprintf(fileID, '%s %9.3f %9.3f %9.3f %9.3f %9.3f  \n', ...
         para_name_Dy(3,:), ... 
       [ mean_lam(3) std_lam(3) para_low(3) para_up(3) inefficiency(3) ]  );  
% fprintf(fileID, '%s %9.3f %9.3f %9.3f %9.3f %9.3f  \n', ...
%        para_name_MS(4,:), ... 
%       [ mean_prob2 std_prob2 para_low(4) para_up(4) inefficiency(4) ]  );    

fprintf(fileID,'\n --------------------------------------------------------------');
fprintf(fileID,'\n --------------------------------------------------------------');
fclose(fileID);


%% ------------------------------------------------------------------------
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
title_name = {'Filtered Posterior lambda (Weight on FA model)'};
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
plot(ti, YY_median(T0_Forecast:Tobs-h_Forecast,1),'LineStyle','-','Color','black', 'LineWidth',2.5); 
title( title_name(1),'FontSize',14 )  
hold off

figure('Name','Predictive Densities of FF, DSGE','Color','w');   
subplot(1,1,1)
    plot(ti', log(PredDen_SW(T0_Forecast:Tobs-h_Forecast)), 'r', 'LineWidth',2)
hold on
    plot(ti', log(PredDen_KK(T0_Forecast:Tobs-h_Forecast)), 'b--', 'LineWidth',2) 
%     plot(ti', log(PredDen_OP(T0_Forecast:Tobs-h_Forecast)),'g-', 'LineWidth',2) 
hold off
        
title( 'Predictive Densities of New Keynesian and Finanical Friction','FontSize',12 );
legend('New Keynesian', 'Financial Friction');

%%



switch parameter
    case 'all'
        save( strcat('./data/DP_WW_mean1_', num2str(horizon),'.mat' ),  'YY_mean', 'YY_median'); 
        save( strcat('./data/save_filter_DP1_', num2str(horizon),'.mat'  ), 'save_filter', 'lambdasim');
    case 'one'
       save( strcat('./data/DP_WW_mean2_', num2str(horizon),'.mat' ),  'YY_mean', 'YY_median'); 
       save( strcat('./data/save_filter_DP2_',  num2str(horizon),'.mat'  ), 'save_filter', 'lambdasim');
end   

est_date = datestr(date);   
switch parameter
    case 'all'
        name = ['./results/Dynamic_Pool_1_',num2str(nsim),'_',est_date,'_', num2str(horizon)];
    case 'one'
       name = ['./results/Dynamic_Pool_2_',num2str(nsim),'_',est_date,'_', num2str(horizon)];
end
   saveas(h_D,name,'fig')
 
 toc