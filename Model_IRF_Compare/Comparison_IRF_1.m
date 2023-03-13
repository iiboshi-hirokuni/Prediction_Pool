
clear all

close all


nirf = 20; 

%------------------------------------------
%  Model with Financial Friction
%------------------------------------------
   Model_FF_post
   yyirf_KK_total = yyirf_total;  % 2014/08/10 C³

%------------------------------------------
%  Model without Financial Friction
%------------------------------------------

    Model_NK_post
    yyirf_DSGE_total = yyirf_total1;  % 2014/08/10 C³

%-------------------------------------------------
%  Drawing Impluse Responses of Structural Shocks
%-------------------------------------------------
    
nshock = 8;  % number of shocks

titlestr = {'Zb  -> ', 'Zg  -> ', 'Zw  -> ', 'Zp  -> ', 'Z \nu  -> ',...
            'Zr  -> ', 'Zz  -> ', 'Z \psi  -> ', 'Z i  -> '};

% ystr = {'Y_t', 'C_t', 'I_t', 'W_t', 'L_t', '\pi_t', 'Pi/P_t',...
%         'R^n_t','E_t \chi_t' ,'B_t','K_t','N_t','Q_t','R^E_t'};

ystr = {'  Output' ,...
        '  Consumption' ,...
			'  Investment' ,...  
			'  Real Wage' ,...
			'  Labor' ,...  
            '  Inflation' ,...
		    '  Investment Price' ,...
             ' Nominal Rate' ,...
             ' Potential GDP',...
             ' Real Borrowing'};

for sh_ind = 1:nshock
     figure('Name','Impluse Responses of Structural Shocks' );
     
     yyirf_KK = squeeze(yyirf_KK_total(:,:,sh_ind));     % 2014/08/10 ’Ç‰Á
     yyirf_DSGE = squeeze(yyirf_DSGE_total(:,:,sh_ind)); % 2014/08/10 ’Ç‰Á
  
       jj=1;
 for j = 1:9    
    if(j<5)||(j==6)||(j==8)
    subplot(2, 3, jj)      
    jj=jj+1;
      hold on
      
      % Impluse Responses with Financial Friction (blue solid line)
      plot(1:nirf, yyirf_KK(:,j),'b-','LineWidth',2)      
      
      % Impluse Responses without Financial Friction (red dashed line)
      plot(1:nirf, yyirf_DSGE(:,j),'r--', 'LineWidth',2)
      
      %  yŽ²‚ª@0‚Æ‚È‚éü
      plot(1:nirf, zeros(nirf,1), 'black' )
      
      title(strcat(titlestr(sh_ind),ystr(j)))
       if j==8
           legend('FA model','NK model');
       end     
      hold off
    end
 end
    
end