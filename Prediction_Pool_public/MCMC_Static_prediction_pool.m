% ------------------------------------------------------------------------
%   Static Predition Pooling 
% ------------------------------------------------------------------------
% clear all;
% clc;

 horizon=2;

%   period = 'bubble'
 period = 'full_sample'
% period = 'Pre_Bubble'
% period = 'Post_Bubble'


varforecast = 6;
k = varforecast+1; 

disp(' Start Staic Prediction Pool')

path('function',path);

 disp('Loading new prediction density');
 path('..\DSGE_KK2\results\FF2011',path);
    load(strcat('FF2011pred_den_', num2str(horizon), '.mat'))
       PredDen_KK = pred_den(:,k);
 path('..\DSGE_SW\results\sw2011',path);     
    load( strcat('sw2011pred_den_', num2str(horizon), '.mat'))
       PredDen_SW = pred_den(:,k);
    BC = csvread('./data/BC_Japan.csv',1,1);   
%     load 'OP2011pred_den.mat'
%        PredDen_OP = pred_den;

       
pd          = 'Normal';
x_old       = 0.0;
lambda_old  = cdf(pd,x_old,0,1);
[Tobs,n]    = size(PredDen_KK); 

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

nsim        = 100000;
nburn       = round(nsim*0.5);

likeold     = -10^(1000);
priornew    = 0;
priorold    = 0;

lambdasim = zeros(nsim,1);

for j = 1:nsim      
    x_new  =  x_old + 0.05*randn(1,1);  
    lambda_new  = x_new;                %lambda_new = cdf(pd,x_new,0,1); 
    
    sum_den     =   0;
    if (x_new > 0) && (x_new < 1)           
        for i = T0_Forecast:Tobs-h_Forecast
            den = (1 - lambda_new)*PredDen_SW(i) + lambda_new*PredDen_KK(i);
            sum_den = sum_den + log(den);
        end   
        likenew = sum_den;
    else
        likenew = -10^(1000);
    end 
    
    % --------------------------------------------------------------------
    %   MH algorithm - compute the acceptance rate 
    % --------------------------------------------------------------------
    postnew = likenew + priornew;      
    postold = likeold + priorold;       
    r = min(1,exp(postnew-postold));   
    if (rand < r)     
        likeold         =  likenew; 
        postold         =  postnew;             
        x_old           =  x_new; 
        lambda_old      =  lambda_new;
        acceptrate(j)   =  1;       
    else
        acceptrate(j)   =   0;
    end            
    lambdasim(j,:) = lambda_old;    
    
     %==     display iteration number ===============     
     itr = num2str(j) ;      
     if mod(j, 10000)==0 
        disp([ itr 'th  iteration.......']);
     end
   %===============================================     
    
    
end        

ti = 1981+(T0_Forecast)/4:0.25:1981+(Tobs-h_Forecast)/4; 
%    
% h_S = figure('Position',[20,20,900,600],'Name','Static Optimal Pools (Weight on FF)','Color','w');
% subplot(1,2,1)
% [density,x1] = ksdensity(lambdasim(nburn:end));
% plot(x1,density,'LineStyle','-','Color','b','LineWidth',2.5);
% title( 'Posterior Density of \lambda','FontSize',12 );
% subplot(1,2,2)
% plot(lambdasim(nburn:end));
% title('Trace of \lambda','FontSize',12 );
% 
%     
% 
% h = figure('Position',[20,20,900,600],'Name','Prediction Score of NK and FF','Color','w');
% hold on
%      area(ti', BC(T0_Forecast:Tobs-h_Forecast)*(-40),'LineStyle','non','FaceColor',[0,0.90,0.90]) 
%     plot(ti', log(PredDen_SW(T0_Forecast:Tobs-h_Forecast)), 'r', 'LineWidth',2)
%     plot(ti', log(PredDen_KK(T0_Forecast:Tobs-h_Forecast)), 'b--', 'LineWidth',2)    
% hold off
%  ylim([-22,-4])
%  ylim       
% title( 'Predictive Densities of KK and SW','FontSize',12 );
% legend('Recessions','NK model', 'FF model');

%--------------------------------------------------------------------------
%   Summary of Posterior Estimates of Lambda  
%--------------------------------------------------------------------------   
mean_lam        = mean(lambdasim(nburn+1:end));
std_lam         = std(lambdasim(nburn+1:end));
sort_para       = zeros(nsim-nburn, 1);
sort_para(:,1)  = sort(lambdasim(nburn+1:end),1);
    
% Calculating of Posterior estimates 
a           =   0.90;   % (1-2*rate) percentage of Credibile interval of Parameters  
rate        =   (1-a)/2;
para_low    =   sort_para(round((nsim-nburn)*rate),1); 
para_up     =   sort_para(round((nsim-nburn)*(1-rate)),1);

para_save   =   lambdasim(nburn:end,1);
p           =   1;
cal_inefficiency    
disp('  mean  s.d. [ lower Band Upper Band ] inefficiency');
disp([ mean_lam std_lam para_low para_up inefficiency ] );

if horizon == 8
     save ./data/mean_lam_8.mat  mean_lam; 
elseif  horizon == 4
     save ./data/mean_lam_4.mat  mean_lam; 
elseif  horizon == 2
     save ./data/mean_lam_2.mat  mean_lam;     
end

%-------------------------------------------------------------
%    output to text file
%-------------------------------------------------------------

para_name_Static = char('\lambda');

est_date = datestr(date);   
result_name = ['./results/Static_pool', est_date , '.txt'];          
fileID = fopen(result_name,'w');
   
fprintf(fileID,'\n---------------------------------------------------------------');
fprintf(fileID,'\n\n                        [ESTIMATION RESULT]');
fprintf(fileID,'\n---------------------------------------------------------------');

fprintf(fileID, '\n Parameter   mean      s.d.     [lower     Upper] inefficiency \n');
fprintf(fileID, '%s %9.3f %9.3f %9.3f %9.3f %9.3f  \n', ...
         para_name_Static(1,:), ...
        [ (mean_lam) (std_lam) (para_low) (para_up) (inefficiency) ] );

fprintf(fileID,'\n --------------------------------------------------------------');
fprintf(fileID,'\n --------------------------------------------------------------');
fclose(fileID);

 est_date = datestr(date);   
% name = ['./results/Static_Pool_',num2str(nsim),'_',est_date];
%         saveas(h_S,name,'fig')

