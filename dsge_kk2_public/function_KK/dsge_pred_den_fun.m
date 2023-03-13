function [forecast] = ...
      dsge_pred_den_fun(para,nforecast,nvar,nshock,t,T1,T0,yy_filter,DD);
%/* INPUT:
%**   para   : parameter vector 
%**   nforecast   : length of forecast
%** OUTPUT:
%**   ssforecast  : impulse response 
%**            the structural shocks.
%**            0 if retcode=0
%**            (nforecast by 25)
%*/

global ZZ YY;


% solve DSGE


    cc00 = 0.25; %0.5;

    impact = cc00* randn(nshock,1);
    
    yyforecast  = zeros(nforecast,nvar);    
    dyyforecast  = zeros(nforecast,nvar);
    
   nstate = size(T0,1);     
    nobs = size(YY,1);
    
     s = yy_filter(t,:)' + [ T0*impact; zeros(size(T0,1),1)]; % 1st period     
%      s =  [ T0*impact; zeros(size(T0,1),1)]; % 1st period     
     ss0  = s(1:nstate,:);        
    
     dyyforecast(1,:) = (ZZ*s)' +DD';     
%     dyyforecast(1,:) = (ZZ*s)' ;

    yyforecast(1,:) = dyyforecast(1,:);
    
    for h = 2:nforecast
        impact = cc00 * randn(nshock,1);
        ss1 = T1*ss0 + T0*impact;    
        s = [ss1;ss0];
        dyyforecast(h,:) = (ZZ*s)'+DD';        
%         yyforecast(t,1:4) = yyforecast(t-1,1:4)+dyyforecast(t,1:4);
%         yyforecast(h,5) = dyyforecast(h,5);
        ss0 = ss1;
    end
    
    forecast = dyyforecast(nforecast,:);
    
       
    yy_est = zeros(nobs,size(YY,2));        
    for i=1:nobs       
        yy_est(i,:) =  (ZZ* yy_filter(i,:)')'+DD'; 
     end    
    
%     yyforecastall(:,1+nvar*(sh_ind-1):nvar*sh_ind) = yyforecast;
    


