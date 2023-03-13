% /*******************************************************/
% /*                                                     */
% /*   Computing Forecast based on                         */
% /*   parameters of DSGE models                           */
% /*                                                     */
% /*******************************************************/

% description: The program converts the DSGE parameters
%              into impulse response functions

disp('start Forecast');

%% forecast period
% nforecast = 20 + 12+20+2;  % From 1989 
% nforecast = 20 + 12;    % From 1991 
nforecast = 20 - 6 ;  % From 1995 

nuse    = 100; % use every nuse observation
nstate = 40;

iparasim = strcat(resupath, runname, '/', runname, 'pa.csv');
fhpara = csvread(iparasim, 1, 0);

 indseq0 = mod(1:size(fhpara,1), nuse);
indseq = (indseq0 ~= 0);
parasim_post  = delif(fhpara,indseq);
nasim = size(parasim_post,1);

hforecast = zeros(nasim,nvar*nforecast);
% hssfilter = zeros(nasim,nobs*nstate);
hssfilter = zeros(nasim,nobs*10);
hrlogh    = zeros(nasim,1);

% iparaest = strcat(postpath,runname,'pm0.csv');
% para = csvread(iparaest);

%*********************************************
%
%  MCMC‚É‚æ‚é DSGE ‚Ì—\‘ª’l
%
%*********************************************

prd_period = nobs-nforecast+1;  

for i = 1:nasim
    % For each element of the block, compute the desired statistic
     para = parasim_post(i,1:npara)';
    
    [ssforecast,yy_est, rlogh] = dsgeforecast_fun(para,nforecast,nvar,nshock,prd_period);
    
    hforecast(i,:) = reshape(ssforecast,nvar*nforecast,1)';
    
%     hssfilter(i,:) = reshape(ssfilter,10*nobs,1)';
    
%     hssfilter(i,:) = reshape(ssfilter,nstate*nobs,1)';
%     hrlogh(i)      = rlogh;

   %==     display iteration number ===============     
     itr = num2str(i) ;      
     nitr = num2str(nasim) ;   
     if mod(i, 100)==0 
        disp([ itr 'th  iteration per ' nitr ' samples' ]);
     end
   %===============================================    

end

forecast_m0 = mean(hforecast, 1);
forecast_s0 = std(hforecast, 1);

hpdprob1_forecast =  0.68; %% 68% interbal
hpdprob2_forecast =  0.95;

drawci1_forecast = []; drawci2_forecast = [];
size_forecast = size(hforecast,2);

for i = 1:size_forecast
    
  drawci1_forecast= [drawci1_forecast, hpdint(hforecast(:,i), hpdprob1_forecast) ];
  
  drawci2_forecast= [drawci2_forecast, hpdint(hforecast(:,i), hpdprob2_forecast) ];
  
end

forecast_hpd_l01 = drawci1_forecast(1,:);
forecast_hpd_h01 = drawci1_forecast(2,:);
forecast_hpd_l02 = drawci2_forecast(1,:);
forecast_hpd_h02 = drawci2_forecast(2,:);

forecast_m     = reshape(forecast_m0,nforecast,nvar);
forecast_hpd_l1 = reshape(forecast_hpd_l01,nforecast,nvar);
forecast_hpd_h1 = reshape(forecast_hpd_h01,nforecast,nvar);
forecast_hpd_l2 = reshape(forecast_hpd_l02,nforecast,nvar);
forecast_hpd_h2 = reshape(forecast_hpd_h02,nforecast,nvar);

%[forecast_hpd_l(:,1),forecast_m(:,1),forecast_hpd_h(:,1)]


ystr = {'  Output' ,...
        '  Consumption' ,...
			' Investment' ,...  
			' Real Wage' ,...
			' Labor' ,...  
            ' Inflation' ,...
		    ' Investment Price' ,...
             ' Nominal Rate' ,...
             ' Real Borrowing',...
             ' Loan rate' ...
             'Net Worth'    };
             

%  h= figure(5500) 
 h= figure('Name', 'forecasting DSGE with FF','Position',[20,20,900,600]) 
 jj=1;
 
  for j = 1:nvar
           
    if(j<5)||(j==6)||(j==8)
    subplot(2, 3, jj)      
    jj=jj+1;
    hold on
    ti2 = 1981+(prd_period-1+1)/4:0.25:1981+(prd_period+nforecast-1)/4;
    
%     hh = area(ti2, [YY(end,j);forecast_hpd_l1(:,j)  0; (forecast_hpd_h1(:,j)-forecast_hpd_l1(:,j))  ] ) ;
     hh = area(ti2, [forecast_hpd_l2(:,j) ...
                    (forecast_hpd_l1(:,j)-forecast_hpd_l2(:,j)) ...
                    (forecast_hpd_h1(:,j)-forecast_hpd_l1(:,j)) ...
                    (forecast_hpd_h2(:,j)-forecast_hpd_h1(:,j))]  ) ;
     
   set(hh(1),'FaceColor',[1 1 1 ])
    set(hh(2),'FaceColor',[1 0.6 0.6])    
    set(hh(4),'FaceColor',[1 0.6 0.6])   
    set(hh(3),'FaceColor',[1 0.4 0.4])  
%     set(hh(2),'FaceColor',[0.753 0.753 0.753])    % Gray
%     set(hh(4),'FaceColor',[0.753 0.753 0.753])   
%     set(hh(3),'FaceColor',[0.553 0.553 0.553])  

    %  Œrü‚ÌÝ’è 
      set(hh,'LineStyle','none');    
    
    
%     plot(ti2, [YY(end,j); forecast_hpd_l1(:,j)] ,'-.r',...
%          ti2, [YY(end,j); forecast_hpd_h1(:,j)] ,'-.r',...
%          ti2, [YY(end,j); forecast_hpd_l2(:,j)] ,'-.r',...
%          ti2, [YY(end,j); forecast_hpd_h2(:,j)] ,'-.r'  )
%      plot(ti2, forecast_m(:,j),'r')

    ti = 1981:0.25:1981+(nobs-1)/4;  
    plot(ti, YY(:,j),'b','LineWidth',2.0);
    
    %     plot(ti2, [YY(end,j); forecast_m(:,j)] ,'b','LineWidth',2.5); 
    plot(ti2, [ forecast_m(:,j)] ,'r','LineWidth',3.5); 

     hold off
     
    title(strcat(ystr(j)))
    if j == 5
        ylim([155 175]);
    end    
    
    end
    
  end
%   w = waitforbuttonpress;

est_date = datestr(date);   
name = ['c:\DSGE\DSGE_KK/results/Forecast_',runname,'_',num2str(nasim),'_',est_date];
        saveas(h,name,'fig')


% figure(5500)  
%   for j = 1:nvar
%     subplot(4, 3, j)
%     
%     ti = 1981:0.25:1981+(nobs-1)/4; 
%         plot(ti, yy_est(:,j),'r-')
%     
%     hold on
%    plot(ti, YY(:,j),'b--')
%      hold off
%      
%     title(strcat(ystr(j)))
%     if j == 1
%         legend('smooth','actual')
%     end    
%   end
