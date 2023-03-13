%====================================================================
% filename:    mcmcdd.g
% description: Computes numerical approximations to the data density
%              based on MCMC output
%====================================================================
% clear all
% clc 

resupath = './results/';
 runname1 = 'FF2011_15000_'
runname = 'FF2011'
runpath = runname;

%densfac = -2150; % level of Log Likelihood; you need to change it, if a model is changed.

hmax    = 20;
n_chains = 1;

thetasim = [];
statsim = [];

for i = 1:n_chains
  %iparasim = strcat(resupath, runpath, '/', runname,'_chain',num2str(i),'_pa.csv');
  iparasim = strcat(resupath, runpath, '/', runname, 'pa.csv');
  pararow  = csvread(iparasim,1,0);
  thetasim = 10*[thetasim;pararow];
  
  %ipostsim = strcat(resupath, runpath, '/', runname,'_chain',num2str(i),'_st.csv');
  ipostsim = strcat(resupath, runpath, '/', runname, 'st.csv');
  statrow  = csvread(ipostsim,1,0);
  statsim = [statsim;statrow];
end

densfac = mean(statsim(:,1)); % level of Log Likelihood;

% --------------------------------
% Method:  modified harmonic mean
% --------------------------------

 p = [0.99 0.95 0.9 0.85 0.8]'; 
% p = [9.9:-0.1:5]'*0.1;
% p = 1-([1:5]'*0.1).^2;

% pcrit = cdfchii(p,ones(size(p, 1),1)*npara);
[npara] = size(thetasim,2);
[ndraws] = size(statsim,1);
pcrit = chi2inv(p,ones(size(p, 1),1)*npara); % statistics toolbox is required

% -------------------------------------------------
% Part 1: Compute posterior moments for parameters

sumtheta  = sum(thetasim,1)';
sumtheta2 = thetasim'*thetasim;

thetahat      = mean(thetasim)'; % sumtheta / ndraws;
% thetasig      = 1.25^(41+2)*cov(thetasim); % ( sumtheta2 / ndraws - thetahat*thetahat');
% thetasig      = diag(100^2.*diag(cov(thetasim)));
thetasig      = cov(thetasim);
thetasiginv   = inv(thetasig);
thetasiglndet = log(det(thetasig));

% -------------------------------------------------
% Part 2:  Compute Density estimate

suminvlike = zeros(size(p, 1),1);
suminvlike2 = zeros(1,1);
laginvlike = zeros(hmax,size(p, 1));
gaminvlike = zeros(hmax,size(p, 1));

invlikemat = [];

for i = 1:ndraws
    lnftheta = -log(p) - 0.5*npara*log(2*pi) - 0.5*thetasiglndet ...
               -0.5*(thetasim(i,:) - thetahat')* ...
               thetasiginv*(thetasim(i,:)' - thetahat);
    
    indtheta = ((thetasim(i,:) - thetahat')* ...
               thetasiginv*(thetasim(i,:)' - thetahat) < pcrit );
    
     invlike  = exp(lnftheta - statsim(i,1) + densfac).*indtheta;
    
    invlikemat = [invlikemat; invlike'];
    
    laginvlike = [invlike'; laginvlike(1:hmax-1,:)];
    for j = 1:hmax
      gaminvlike(j,:) = gaminvlike(j,:) + laginvlike(j,:).*invlike';
    end    
    
    suminvlike = suminvlike + invlike;
    
%    Marginal likelihood (Newton and Raftery,1994)
	 invlike2 = exp(- statsim(i,2)+densfac); 
     suminvlike2 = suminvlike2 + invlike2;
end

meaninvlike = suminvlike/ndraws;

meaninvlike2 = suminvlike2/ndraws;

for j = 1:hmax
  gaminvlike(j,:) = gaminvlike(j,:)/ndraws - (meaninvlike.^2)';
end

% -------------------------------------------------
% Compute Standard Errors

suminvlikeerror0 = gaminvlike(1,:);

j = 2;
while j <= hmax
   suminvlikeerror0 = suminvlikeerror0 + 2*gaminvlike(j,:)*(1-(j-1)/hmax);
   j = j+1;
end
serror = 100*sqrt(suminvlikeerror0/ndraws)./meaninvlike';

pym3 = densfac-log(meaninvlike);   % log marginal likelihood
pym4 = densfac-log(meaninvlike2);   % log marginal likelihood (Newton and Raftery,1994)

disp('Marginal Likelihood (Geweke 1999)')
disp(sprintf('%20s%10.5g%10.5g%10.5g%10.5g%10.5g',...
             'pcrit %',p(1),p(2),p(3),p(4),p(5)))
disp(sprintf('%20s%10.5g%10.5g%10.5g%10.5g%10.5g',...
             'Log Marginal lkh',...
             pym3(1),pym3(2),pym3(3),pym3(4),pym3(5)))
disp(sprintf('%20s%10.5g%10.5g%10.5g%10.5g%10.5g',...
             'Simulation Error(%)',...
             serror(1),serror(2),serror(3),serror(4),serror(5)))
disp(sprintf('%20s%10.5g%', 'Average Across p =', mean(pym3)))
disp(sprintf('%20s%10.5g%', 'Marginal likelihood (Newton and Raftery, 1994) =', pym4) );

est_date = datestr(date);   
result_name = ['c:\DSGE\DSGE_KK/results/Marginal_Likelihood_',runname,'_',num2str(ndraws), '_', est_date , '.txt'];          
fileID = fopen(result_name,'w');
   
fprintf(fileID,'\n\n                        [ESTIMATION RESULT]');
fprintf(fileID,'\n-------------------------------------------------------------------------------');
fprintf(fileID,'\n----------------------------------------------------------------------------- \n');
fprintf(fileID,'\n %20s%10.5g%10.5g%10.5g%10.5g%10.5g',...
             'pcrit %',p(1),p(2),p(3),p(4),p(5));
fprintf(fileID,'\n %20s%10.5g%10.5g%10.5g%10.5g%10.5g',...
             'Log Marginal lkh ',...
             pym3(1),pym3(2),pym3(3),pym3(4),pym3(5));
fprintf(fileID,'\n %20s%10.5g%10.5g%10.5g%10.5g%10.5g',...
             'Simulation Error(%)',...
             serror(1),serror(2),serror(3),serror(4),serror(5));

fprintf(fileID,'\n-------------------------------------------------------------------------------');
fprintf(fileID,'\n-------------------------------------------------------------------------------');              
fprintf(fileID,'\n %20s%10.5g%', 'Average Across p =', mean(pym3) );
fprintf(fileID,'\n-------------------------------------------------------------------------------');
fprintf(fileID,'\n-------------------------------------------------------------------------------');
fprintf(fileID,'\n-------------------------------------------------------------------------------');
fprintf(fileID,'\n-------------------------------------------------------------------------------');              
fprintf(fileID,'\n %20s%10.5g%', ' Marginal likelihood (Newton and Raftery) =', mean(pym4) );
fprintf(fileID,'\n-------------------------------------------------------------------------------');
fprintf(fileID,'\n-------------------------------------------------------------------------------');
fclose(fileID);

%"Marginal Data Density"; 
%
%pcrit_names = 10 | 20| 30 | 40 | 50 |60|70| 80| 90   ;
%
%pcrit_names~pym3;
%
%"";
%"Average Across p     " meanc(pym3);
%
%"";
%"Simulation Error [Percent]" suminvlikeerror';
%
%title("Data Density and Marginal Likelihood");
%xy(pcrit_names, pym3~meanc(pym3)*ones(9,1));

