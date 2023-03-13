%******************************************************/
%                                                     */
% Computing Moments based on parameters of DSGE model */
%                                                     */
%******************************************************/

% filename:    nkbc_ac.g
% description: The program converts the DSGE parameters
%              into unconditional moments

nac = nvar + nvar*(nvar+1)/2 + nvar^2;
nvd = nvar * 9;

nn = 4;   % numbers of Horizon

% Open files with Parameter Draws

iparasim = strcat(resupath, runname, '/', runname, 'pa.csv');
fhpara = csvread(iparasim, 1, 0);

% Define some Constants

nuse    = 1;      % use every nuse observation

hacdraw = [];
hvddraw = [];

indseq0 = mod(1:nsim, nuse);
indseq = (indseq0 ~= 0);
parasim = delif(fhpara,indseq);
nasim = size(parasim,1);

for j = 1:nasim
    
        % For each element of the block, compute the desired statistic
        para = parasim(j,1:npara)';
        
        [GAM0_c, GAM0_i, GAM0_q, GAM0_l, GAM0_w, GAM0_z, GAM0_p, GAM0_g, GAM0_m, ...
          GAM1_c, GAM1_i, GAM1_q, GAM1_l, GAM1_w, GAM1_z, GAM1_p, GAM1_g, GAM1_m ...
            ] = nkbcmom(para, nn);
        
        GAM0 = GAM0_c + GAM0_i + GAM0_q + GAM0_l ...
            + GAM0_w + GAM0_z + GAM0_p + GAM0_g + GAM0_m;
        GAM1 = GAM1_c + GAM1_i + GAM1_q + GAM1_l ...
            + GAM1_w + GAM1_z + GAM1_p + GAM1_g + GAM1_m;
        
        % Variance Decompositions
        vd  = [(diag(GAM0_c)./diag(GAM0)), (diag(GAM0_i)./diag(GAM0)),...
               (diag(GAM0_q)./diag(GAM0)), (diag(GAM0_l)./diag(GAM0)),...
               (diag(GAM0_w)./diag(GAM0)), (diag(GAM0_z)./diag(GAM0)),...
               (diag(GAM0_p)./diag(GAM0)), (diag(GAM0_g)./diag(GAM0)),...
               (diag(GAM0_m)./diag(GAM0))];
        %vdsim[j,.] =  vec(vd')';
        hvddraw(j,:) = reshape(vd', nvd, 1);
        
        % Autocovariances
        corr0  = corrvc(GAM0);
        stddev = sqrt(diag(GAM0));
        corr1  = GAM1 ./ (stddev*stddev');
        
        hacdraw(j,:) = [stddev' , vech(corr0)' ,reshape(corr1',nvar*nvar, 1)'];

end

oac  = strcat(resupath, runname, '/', runname, 'ac.csv');
csvwrite(oac, hacdraw)

ovd  = strcat(resupath, runname, '/', runname, 'vd.csv');
csvwrite(ovd, hvddraw)

drawmean = mean(hacdraw, 1);
drawstdd = std(hacdraw, 1);
drawcov = cov(hacdraw);

hpdprob = 0.90;
drawci = [];

for i = 1:nvd
  drawci= [drawci, hpdint(hacdraw(:,i), hpdprob)];
end

draw_names = { 'Y_c','Y_i','Y_q','Y_L','Y_w','Y_z','Y_p','Y_g','Y_m',...
               'Pi_c','Pi_i','Pi_q','Pi_L','Pi_w','Pi_z','Pi_P','Pi_g','Pi_m',...
               'W_c','W_i','W_q','W_L','W_w','W_z','W_p','W_g','W_m',...
               'I_c','I_i','I_q','I_L','I_w','I_z','I_p','I_g','I_m',...
               'C_c','C_i','C_q','C_L','C_w','C_z','C_p','C_g','C_m',...
               'R_c','R_i','R_q','R_L','R_w','R_z','R_p','R_g','R_m',...
               'L_c','L_i','L_q','L_L','L_w','L_z','L_p','L_g','L_m'};

disp('');
disp('===========================================');

disp(sprintf('%15s%15s%15s%15s%15s',...
        'Parameter','Mean','StdD','CI(Low)','CI(High)'));
for i = 1:length(draw_names)
    disp(sprintf('%15s%15.5g%15.5g%15.5g%15.5g',...
       char(draw_names(i)),drawmean(i),drawstdd(i),drawci(1,i),drawci(2,i)));
end

