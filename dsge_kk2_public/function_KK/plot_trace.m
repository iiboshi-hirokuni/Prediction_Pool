% 
% 
%  Plot traces of sampling parameters
% 
% 
% 



para_names_p = {'\sigma','\theta','\chi','inv_{\zeta}',...
     '\mu','\phi_{o/y}','\gamma_w','\xi_w', ...
     '\gamma_p','\xi_p','\phi_r','\phi_{\pi}','\phi_y','{z}^{*}','\psi_{bar}', ...
     '\eta','\nu_k','\mu_E','r^E_{bar}', ...
     '\rho_b','\rho_g','\rho_w','\rho_p','\rho_r','\rho_{\nu}','\rho_z','\rho_{\psi}', ...
     '\rho_{efp}','\rho_{nw}', ...
     '\sigma_b','\sigma_g','\sigma_w','\sigma_p','\sigma_r', ...
     '\sigma_{nu}','\sigma_z','\sigma_{\psi}','\sigma_{efp}','\sigma_{nw}'};
 


%% Output posterior file
  nuse    = 1;
  iparasim = strcat(resupath, runname, '/', runname, 'pa.csv');
   fhpara = csvread(iparasim, 1, 0);
      
      indseq0 = mod(1:(nblock-nburnin_block)*nsim, nuse);
      indseq = (indseq0 ~= 0);
      parasim_post  = delif(fhpara,indseq);
      
   nasim = size(parasim_post,1);




%% draw graph both of Prior and posterior densities
    k = 0;
    
 for i =1:1:npara
    
  if mod(i,3*4)==1
     k = k +1 ;
     figure(2500+k*10)     
  end
 
 subplot(3,4,i-(k-1)*3*4) % 3(çs)Å~4(óÒ)
 
%      [density0,x0]  = ksdensity(parasim_pri(1:end,i));
%      [density,x]  = ksdensity(parasim_post(1:end,i));
     
           plot(parasim_post(:,i)','LineStyle','-','Color','b',...
          'LineWidth',1.0);
      
%      if both_plot == 1 
%       hold on
%           plot(x0,density0,'LineStyle','--','Color','r',...
%          'LineWidth',2.5);
%        hold off 
%      end
      
     title( char(para_names_p(i)),'FontSize',12) ; 
  
 end





