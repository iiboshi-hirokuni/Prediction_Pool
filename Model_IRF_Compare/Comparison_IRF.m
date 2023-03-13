close all

nirf = 20; 

%------------------------------------------
%  Model with Financial Friction
%------------------------------------------
   Model_FF_post
   yyirf_KK_total = yyirf_total;  % 2014/08/10 èCê≥

%------------------------------------------
%  Model without Financial Friction
%------------------------------------------

    Model_NK_post
    yyirf_DSGE_total = yyirf_total1;  % 2014/08/10 èCê≥

%-------------------------------------------------
%  Drawing Impluse Responses of Structural Shocks
%-------------------------------------------------
    
nshock = 9;  % number of shocks

titlestr = {'Zb  -> ', 'Zg  -> ', 'Zw  -> ', 'Zp  -> ', 'Z \nu  -> ',...
            'Zr  -> ', 'Zz  -> ', 'Z \psi  -> ', 'Z i  -> '};

ystr = {'Y_t', 'C_t', 'I_t', 'W_t', 'L_t', '\pi_t', 'Pi/P_t',...
        'R^n_t','E_t \chi_t' ,'B_t','K_t','N_t','Q_t','R^E_t'};


for sh_ind = 1:nshock
     figure('Position',[20,20,900,600],'Name','Impluse Responses of Structural Shocks' );
     
     yyirf_KK = squeeze(yyirf_KK_total(:,:,sh_ind));     % 2014/08/10 í«â¡
     yyirf_DSGE = squeeze(yyirf_DSGE_total(:,:,sh_ind)); % 2014/08/10 í«â¡
     
 for j = 1:9
      subplot(3, 3, j)
      hold on
      
      % Impluse Responses with Financial Friction (blue solid line)
      plot(1:nirf, yyirf_KK(:,j),'b-','LineWidth',2)      
      
      % Impluse Responses without Financial Friction (red dashed line)
      plot(1:nirf, yyirf_DSGE(:,j),'r--', 'LineWidth',2)
      
      %  yé≤Ç™Å@0Ç∆Ç»ÇÈê¸
      plot(1:nirf, zeros(nirf,1), 'black' )
      
      title(strcat(titlestr(sh_ind),ystr(j)))
       if j==9
           legend('FA Model','NK model');
       end     
      hold off
 end
    
end