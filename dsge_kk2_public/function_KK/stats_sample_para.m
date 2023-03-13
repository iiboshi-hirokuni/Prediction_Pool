



para_names_p = char('\sigma','\theta','\chi','inv_{\zeta}',...
     '\mu','\phi_{o/y}','\gamma_w','\xi_w', ...
     '\gamma_p','\xi_p','\phi_r','\phi_{\pi}','\phi_y','{z}^{*}','\psi_{bar}', ...
     '\eta','\nu_k','\mu_E','r^E_{bar}', ...
     '\rho_b','\rho_g','\rho_w','\rho_p','\rho_r','\rho_{\nu}','\rho_z','\rho_{\psi}', '\rho_{i}',...
     '\rho_{efp}','\rho_{nw}', ...
     '\sigma_b','\sigma_g','\sigma_w','\sigma_p','\sigma_r', ...
     '\sigma_{nu}','\sigma_z','\sigma_{\psi}','\sigma_{i}','\sigma_{efp}','\sigma_{nw}');
 


%% Output posterior file
  nuse    = 10;
  iparasim = strcat(resupath, runname, '/', runname, 'pa.csv');
   fhpara = csvread(iparasim, 1, 0);

%   indseq0 = mod(1:(nblock-nburnin_block)*nsim/nthin, nuse);
    indseq0 = mod(1:size(fhpara,1), nuse);
      indseq = (indseq0 ~= 0);
      parasim_post  = delif(fhpara,indseq);
   
   nasim = size(parasim_post,1);
   npara = size(parasim_post,2);
   
   % Calculating of Posterior estimates 

a  =0.90;   % (1-2*rate) percentage of Credibile interval of Parameters  
rate = (1-a)/2;


% calculation of posterior estimates of parameters
%

sort_para=zeros( nasim, npara );


  for i=1:1:npara
     sort_para(:,i) = sort(parasim_post(:,i),1);
  end


%  The Fisrt equation
  para_low = sort_para(round(( nasim )*rate),:); 
  para_up  = sort_para(round(( nasim )*(1-rate)),:);

iBm = min([500, nasim/2]); 
 
est_date = datestr(date);   
result_name = ['c:\DSGE\DSGE_KK/results/ESTIMATE_para_',runname,'_',num2str(nasim), '_', est_date , '.txt'];          
fileID = fopen(result_name,'w');
   
fprintf(fileID,'\n\n                        [ESTIMATION RESULT]');
fprintf(fileID,'\n----------------------------------');
fprintf(fileID,'------------------------------------');
fprintf(fileID,'\nParameter         Mean        Stdev     ');
fprintf(fileID,'95%%Low     95%%Up    Geweke     Inef.');
fprintf(fileID,'\n----------------------------------');
fprintf(fileID,'--------------------------------------------\n');

for i = 1:npara;
% fprintf( '%s ',  para_names_p(i)  );
fprintf(fileID,'%s %10.4f  %10.4f %9.3f %9.3f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim_post(:,i)) std(parasim_post(:,i)) para_low(i) para_up(i) ...
      fGeweke(parasim_post(:,i), iBm), ...
      ftsvar(parasim_post(:,i), iBm)/var(parasim_post(:,i)) ]  );
end


fprintf(fileID,'-----------------------------------');
fprintf(fileID,'-----------------------------------');
fclose(fileID);
