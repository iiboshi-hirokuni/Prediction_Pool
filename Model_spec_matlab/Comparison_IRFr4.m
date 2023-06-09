close all
clear all   %2014.8.15中村修正
nirf = 20; %2014.8.15中村修正

%------------------------------------------
%  Model with Financial Friction
%------------------------------------------
   Modelspec_KKr3
   yyirf_KK_total = yyirf_total;  % 2014/08/10 修正

%------------------------------------------
%  Model without Financial Friction
%------------------------------------------

%    Modelspec_DSGE
    PLANE_VANILAr3   % <---  ここを入れ替えてください
    yyirf_DSGE_total = yyirf_total1;  % 2014/08/10 修正
    
%------------------------------------------
%  Small Open Model without Financial Friction
%------------------------------------------

%    Modelspec_DSGE
    Small_Open_KKr1   % <---  ここを入れ替えてください
    yyirf_SM_total = yyirf_total2;  % 2014/08/10 修正    
    
%------------------------------------------
%  Unemployment Rate (Gali et al. 2011 )
%------------------------------------------

%    Modelspec_DSGE
    Gali_Unemp   % <---  ここを入れ替えてください
    yyirf_UR_total = yyirf_total3;  % 2014/08/10 修正   
    

%-------------------------------------------------
%  Drawing Impluse Responses of Structural Shocks
%-------------------------------------------------
    
nshock = 8;  % number of shocks


for sh_ind = 1:nshock
     figure('Name','Impluse Responses of Structural Shocks' );
     
     yyirf_KK = squeeze(yyirf_KK_total(:,:,sh_ind));     % 2014/08/10 追加
     yyirf_DSGE = squeeze(yyirf_DSGE_total(:,:,sh_ind)); 
     yyirf_SM = squeeze(yyirf_SM_total(:,:,sh_ind));  
     yyirf_UR = squeeze(yyirf_UR_total(:,:,sh_ind)); 
     
 for j = 1:16 %2014.8.15中村修正
      subplot(4, 4, j) %2014.8.15中村修正
      hold on
      
      if j <= 14
      % Impluse Responses with Financial Friction (blue solid line)
      plot(1:nirf, yyirf_KK(:,j),'b-','LineWidth',2)      
      
      % Impluse Responses without Financial Friction (red dashed line)
      plot(1:nirf, yyirf_DSGE(:,j),'r--', 'LineWidth',2)
      end
      
      if j <= 15
      % Impluse Responses of Small Open (green dashed line)
      plot(1:nirf, yyirf_SM(:,j),'g--', 'LineWidth',2)   
      end
      
      if j <= 14 | j ==16
      % Impluse Responses of Gali (black dashed line)
      plot(1:nirf, yyirf_UR(:,j),'black--', 'LineWidth',1.0) 
      end
      
      %  y軸が　0となる線
      plot(1:nirf, zeros(nirf,1), 'black' )
      
      title(strcat(titlestr(sh_ind),ystr(j)))
       if j==1
           legend('Financial Friction','DSGE','Small Open', 'Unemp');
       end     
      hold off
 end
    
end