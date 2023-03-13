% /*******************************************************/
% /*                                                     */
% /*    Computing Prediction Pool based on               */
% /*    parameters of DSGE models                        */
% /*                                                     */
% /*******************************************************/



horizon = 8;


% target of forecasting
j1 = 1;     % output=1,
j2 = 2;     % cons = 2
j3 = 3;     % inv
j4 = 4;     % real wage
j5 = 6;     %  inflation=6
j6 = 8;     % Nominal Interest Rate= 8

var5 = [1 2 3 4 6 8];

nuse    = 10; % use every nuse observation
nstate = 40;

% iparasim = strcat(resupath, runname, '/', runname, 'pa.csv');
iparasim = strcat(resupath, runname, '/', runname, '_15000_pa.csv');
fhpara = csvread(iparasim, 1, 0);

indseq0 = mod(1:size(fhpara,1), nuse);
indseq = (indseq0 ~= 0);
parasim_post  = delif(fhpara,indseq);
nasim = size(parasim_post,1);


%*********************************************
%
%  DSGE の prediction density の計算
%
%*********************************************

varforcast = 6;

pred_den = zeros(nobs, varforcast+2);
den = zeros(varforcast+2,1);

save_diff = zeros(nobs, varforcast);
diff = zeros(varforcast,1);


 for i = 1:nasim
      
     para = parasim_post(i,1:npara)';
     
     [T1,TC,T0,RC] = dsgesolv(para);
    
     [rloglh,yy_filter,DD] = Kalman_dsge(para,YY,nshock,ZZ,T1,TC,T0);

  for t = horizon+1:nobs
    
     for h = 1:horizon      
    
    [yy_f] = dsge_pred_den_fun(para,h,nvar,nshock,t-h,T1,T0,yy_filter,DD);
    
    sigma = eye(varforcast);
    sigma(1,1) = 1/gamrnd(1, 3.0^2);  % output
    sigma(2,2) = 1/gamrnd(1, 3.0^2);  % cons
    sigma(3,3) = 1/gamrnd(1, 3.0^2);  % inv
    sigma(4,4) = 1/gamrnd(1, 3.0^2);  % wage
    sigma(5,5) = 1/gamrnd(1, 3.0^2);  % inflation
     sigma(6,6) = 1/gamrnd(1, 3.0^2);  % inflation
    
%     den(varforcast+1) = mvnpdf( [ YY(t,j1) YY(t,j2) YY(t,j3) YY(t,j4) YY(t,j5) ],...
%                    [yy_f(j1) yy_f(j2) yy_f(j3) yy_f(j4) yy_f(j5)  ] , sigma );
               
   den(varforcast+1) = mvnpdf( [ YY(t,j1) YY(t,j2) YY(t,j3) YY(t,j4) YY(t,j5) YY(t,j6) ],...
                   [yy_f(j1) yy_f(j2) yy_f(j3) yy_f(j4) yy_f(j5) yy_f(j6) ] , sigma );              
    
    pred_den(t,varforcast+1) = pred_den(t,varforcast+1)+ den(varforcast+1);  %% horizon の合計
  
    
    %% forecasting without investment
    den(varforcast+2) = mvnpdf( [ YY(t,j1) YY(t,j2) YY(t,j4) YY(t,j5) YY(t,j6) ],...
                   [yy_f(j1) yy_f(j2) yy_f(j4) yy_f(j5) yy_f(j6) ] , ...
                    diag([sigma(1,1) sigma(2,2) sigma(4,4) sigma(5,5) sigma(6,6)]) );              
    
    pred_den(t,varforcast+2) = pred_den(t,varforcast+2)+ den(varforcast+2);  %% horizon の合計
    
      for j = 1:varforcast
         den(j) = mvnpdf( YY(t,var5(j)), yy_f(var5(j)), sigma(j,j) ); 
         pred_den(t,j) = pred_den(t,j) + den(j);  %% horizon の合計
         
%          if den(j) > mvnpdf(YY(t,var5(j)),YY(t,var5(j))+2*sqrt(sigma(j,j)) )
            diff(j) =  yy_f(var5(j)) - YY(t,var5(j)); %% forecast - actual
%             diff(j) =  yy_f(var5(j)) ;
%            if (j==6)&&(t>45)              
%                disp([ num2str(1981+0.25*t) ', ' num2str(yy_f(var5(j))) ', '...
%                    num2str(YY(t,var5(j))) ', ' num2str(yy_f(var5(j))-YY(t,var5(j))) ]); 
%            end    
             save_diff(t,j) = save_diff(t,j) + diff(j);
%          end
         
      end
    
    end

  end
  
  if mod(i/nasim, 0.1)==0
     disp(sprintf('complete parcentage of nsim:    %f\n',i/nasim ));
  end
     
end

  pred_den = pred_den./(nasim*horizon);
  
 figure(7000) 
 ti = 1981:0.25:1981+(nobs-1)/4;  
    plot(ti,pred_den,'LineWidth',3);
    legend('y','cons','inv','wage','inf','Nominal rate','Total');

opred = strcat(resupath, runname, '/', runname, 'pred_den_',num2str(horizon),'.mat'); 
save(opred, 'pred_den'); 

save_diff = save_diff/nasim/horizon;
% save_diff = save_diff/nasim;

opred = strcat(resupath, runname, '/', runname, 'pred_diff.mat'); 
 save(opred, 'save_diff'); 

var_title ={'output','cons','inv','wage','inf','Nominal rate' };

% figure(8000)
% for i = 1:6
%   subplot(6,1,i)
%   hold on
%     plot(ti,save_diff(:,i)/nasim/horizon )
%     plot(ti,zeros(size(ti,2),1),'r--' )
%   hold off 
%    title(var_title(i))
% end

