% /*******************************************************/
% /*                                                     */
% /* Computing Historical Decomposition based on                         */
% /* parameters of DSGE models                           */
% /*                                                     */
% /*******************************************************/

% description: The program converts the DSGE parameters
%              into impulse response functions

nbos        = size(YY,1);
% nshock      = 11;
nhistdecomp = nbos*nshock; 
nuse    = 20; % use every nuse observation
nstate = 40;

iparasim = strcat(resupath, runname, '/', runname, 'pa.csv');
fhpara = csvread(iparasim, 1, 0);

indseq0 = mod(1:(nblock-nburnin_block)*nsim, nuse);
indseq = (indseq0 ~= 0);
parasim_post  = delif(fhpara,indseq);
nasim = size(parasim_post,1);

hhistdecomp = zeros(nasim,nvar*nhistdecomp);
% hssfilter = zeros(nasim,nobs*nstate);
hssfilter = zeros(nasim,nobs*10);
hrlogh    = zeros(nasim,1);

leverage_save = zeros(nobs,nasim);

% iparaest = strcat(postpath,runname,'pm0.csv');
% para = csvread(iparaest);

for i = 1:nasim
    % For each element of the block, compute the desired statistic
     para = parasim_post(i,1:npara)';
    
    [yyhistdecomp,yy_smooth1,yy_smooth2,DD,leverage] =...
        dsgehistdecomp_fun(para,nhistdecomp,nvar,nshock);
    
    hhistdecomp(i,:) = reshape(yyhistdecomp,nvar*nhistdecomp,1)';
    
    leverage_save(:,i) =  leverage;

end

histdecomp_m0 = mean(hhistdecomp, 1);
histdecomp_s0 = std(hhistdecomp, 1);

hpdprob1_histdecomp =  0.68; %% 68% interbal
hpdprob2_histdecomp =  0.95;

drawci1_histdecomp = []; drawci2_histdecomp = [];
size_histdecomp = size(hhistdecomp,2);

for i = 1:size_histdecomp
    
  drawci1_histdecomp= [drawci1_histdecomp, hpdint(hhistdecomp(:,i), hpdprob1_histdecomp) ];
  
  drawci2_histdecomp= [drawci2_histdecomp, hpdint(hhistdecomp(:,i), hpdprob2_histdecomp) ];
  
end

histdecomp_hpd_l01 = drawci1_histdecomp(1,:);
histdecomp_hpd_h01 = drawci1_histdecomp(2,:);
histdecomp_hpd_l02 = drawci2_histdecomp(1,:);
histdecomp_hpd_h02 = drawci2_histdecomp(2,:);

histdecomp_m     = reshape(histdecomp_m0,nobs,nvar*nshock);
histdecomp_hpd_l1 = reshape(histdecomp_hpd_l01,nobs,nvar*nshock);
histdecomp_hpd_h1 = reshape(histdecomp_hpd_h01,nobs,nvar*nshock);
histdecomp_hpd_l2 = reshape(histdecomp_hpd_l02,nobs,nvar*nshock);
histdecomp_hpd_h2 = reshape(histdecomp_hpd_h02,nobs,nvar*nshock);

%[histdecomp_hpd_l(:,1),histdecomp_m(:,1),histdecomp_hpd_h(:,1)]


ystr = {'  Output' , '  Consumption' , 'Investment' , ' Real Wage' ,...
			' Labor' , ' Inflation' , ' Investment Price' , ' Nominal Rate' ,...
             ' Potential GDP', ' Real Borrowing'};
         
%          Epssilon_t = [epssilon_b; epssilon_g; epssilon_w; epssilon_p; epssilon_nu; ...
%              epssilon_r; epssilon_z; epssilon_pssi; epssilon_i; epssilon_efp; epssilon_nw];
         
 titlestr = {'Preference Shock: Z^b ' ,...
           'Exg-Goods Demand Shock: Z^g' ,...
		   'Wage Markup Shock: Z^w ' ,...
		   'Price Markup Shock: Z^p ' ,...
		   'Investment Price Shock: Z^{\nu}' ,...
		   'Monetary Policy Shock: Z^r' ,...
		   'Productivity Shock: Z^z ' ,... 
		   'Investment Shock: Z^{\psi}' ,...
           'Investment Shock: Z^{i}' ,...
		   'Financial Shock: Z^{efp}' ,...
		   'Net Worth Shock: Z^{nw}'};  
  
       
 ti = 1981.5:0.25:1981.5+(nobs-2)/4;    
 
 for i = 1:nvar
   figure(6000+10*i)    
    
    subplot(4, 3, 1)
     hist_total = zeros(nobs,1); 
    for k = 1:nshock
       hist_total = hist_total + histdecomp_m(:,i+(k-1)*nvar);
    end
          plot(ti,hist_total(2:end),'b');
       hold on
           plot(ti,YY(2:end,i),'r--');
       hold off
        title( [ strcat(ystr(i)) ] )
         legend( 'estimated', 'actual' )
         
    subplot(4, 3, 12) 
%          plot(ti,YY(2:end,i)-hist_total(2:end),'b-');
          area(ti,YY(2:end,i)-hist_total(2:end));
     hold on
        plot(ti,YY(2:end,i),'r-','LineWidth', 1.5);
     hold off
         
         title( 'Impact of Initial value'  )
         legend( 'initial' )
    
    
    
  for j = 1:nshock
    subplot(4, 3, j+1)
    
   
%      plot(ti,histdecomp_m(2:end,i+(j-1)*nvar),'b')
       area(ti,histdecomp_m(2:end,i+(j-1)*nvar))
    
    hold on
         plot(ti,YY(2:end,i),'r-','LineWidth', 1.5);
%       plot(ti, histdecomp_hpd_l1(2:end,i+(j-1)*nshock),'-.r')
%       plot(ti, histdecomp_hpd_h1(2:end,i+(j-1)*nshock),'-.r' )
%       plot(ti, zeros(size(ti,2),1),'k-' )
     hold off
     
    title( [ strcat(titlestr(j)) ] )
    if j == nshock
     legend( 'contribution', 'observed series' )
    end
  end  
  
  
  end
%   w = waitforbuttonpress;


figure(5500)  
  for j = 1:nvar
    subplot(4, 3, j)
    
    ti = 1981.5:0.25:1981+(nobs-1)/4; 
        plot(ti, yy_smooth1(3:end,j),'r-')
    
    hold on
     plot(ti, yy_smooth2(3:end,j),'k:')
     plot(ti, YY(3:end,j),'b--')  
     hold off
     
    title(strcat(ystr(j)))
    if j == 1
        legend('state-smooth','dis-smooth','actual')
    end    
  end
 
  figure(5600)  
  
    subplot(1, 1, 1)
    
    ti = 1981.5:0.25:1981+(nobs-1)/4; 
    
%     hold on
     plot(ti,  mean(leverage_save(3:end,:),2),'k-')
%     hold off

    

