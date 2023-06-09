close all

nirf = 20; 

%------------------------------------------
%  Model with Financial Friction
%------------------------------------------
   Modelspec_KK
   yyirf_KK_total = yyirf_total;  % 2014/08/10 修正

%------------------------------------------
%  Model without Financial Friction
%------------------------------------------

%    Modelspec_DSGE
    PLANE_VANILAr2   % <---  ここを入れ替えてください
    yyirf_DSGE_total = yyirf_total1;  % 2014/08/10 修正

%-------------------------------------------------
%  Drawing Impluse Responses of Structural Shocks
%-------------------------------------------------
    
nshock = 8;  % number of shocks


for sh_ind = 1:nshock
     figure('Name','Impluse Responses of Structural Shocks' );
     
     yyirf_KK = squeeze(yyirf_KK_total(:,:,sh_ind));     % 2014/08/10 追加
     yyirf_DSGE = squeeze(yyirf_DSGE_total(:,:,sh_ind)); % 2014/08/10 追加
     
 for j = 1:9
      subplot(3, 3, j)
      hold on
      
      % Impluse Responses with Financial Friction (blue solid line)
      plot(1:nirf, yyirf_KK(:,j),'b-','LineWidth',2)      
      
      % Impluse Responses without Financial Friction (red dashed line)
      plot(1:nirf, yyirf_DSGE(:,j),'r--', 'LineWidth',2)
      
      %  y軸が　0となる線
      plot(1:nirf, zeros(nirf,1), 'black' )
      
      title(strcat(titlestr(sh_ind),ystr(j)))
       if j==1
           legend('Financial Friction','DSGE');
       end     
      hold off
 end
    
end