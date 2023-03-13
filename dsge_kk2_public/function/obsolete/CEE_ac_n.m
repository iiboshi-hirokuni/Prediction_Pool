%******************************************************/
%                                                     */
% Computing Moments based on parameters of DSGE model */
%                                                     */
%******************************************************/

% filename:    nkbc_ac.g
% description: The program converts the DSGE parameters
%              into unconditional moments

%  __output = 1;
%  _fcmptol = 1E-10;

%nvar = 7;  % num of  observed variables
nac = nvar + nvar*(nvar+1)/2 + nvar^2;
nvd = nvar * 9;

nn = 4;   % numbers of Horizon

% Open files with Parameter Draws

iparasim = strcat(resupath, runname, '/', runname, 'pa.csv');
fhpara = csvread(iparasim, 1, 0);

% Define some Constants

%nblock  = 1000;   % block size
nuse    = 1;      % use every nuse observation
%nsim

% Initialize Output files
%create fhac = ^oac with AC, nac, 8;
%create fhvd = ^ovd with VD, nvd, 8;
 
%/* Initialization
%*/ 
%eofloop = 0;
%loopct  = 1;
%
%do until eofloop; 

hacdraw = [];
hvddraw = [];

for i = 1:nblock

    % extract sampled parameters
    parasim = fhpara((1:nsim)+(i-1)*nsim,:);
    
    indseq0 = mod(1:nsim, nuse);
    indseq = (indseq0 ~= 0);
    
    parasim = delif(parasim,indseq);
    nasim = size(parasim,1);
    
    acsim = zeros(nasim,nac);
    vdsim = zeros(nasim,nvd);
    
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
        vdsim(j,:) = reshape(vd', nvd, 1);
        
        % Autocovariances
        corr0  = corrvc(GAM0);
        stddev = sqrt(diag(GAM0));
        corr1  = GAM1 ./ (stddev*stddev');
        
        acsim(j,:) = [stddev' , vech(corr0)' ,reshape(corr1',nvar*nvar, 1)'];
        
    end
    hacdraw = [hacdraw; acsim];
    hvddraw = [hvddraw; vdsim];
    
end

oac  = strcat(resupath, runname, '/', runname, 'ac.csv');
csvwrite(oac, hacdraw)

ovd  = strcat(resupath, runname, '/', runname, 'vd.csv');
csvwrite(ovd, hvddraw)

% oldfilename:    vdpmom.g
% description: Compute posterior mean and std and CI based on
%              output of posterior simulator
%

%/******************************************************************
%**         Load Parameter Draws from MH Output
%*/
%#include c:\dge_sw\pathspec.g;
%#include c:\dge_sw\ceespec.g;
%
%mhrun  = "2";    @  2: SW model, 3: Case A, 4: Case B, 5: Case C   @
%
%lpath  = resupath $+ "\\mhrun" $+ mhrun;
%opath  = resupath $+ "\\mhrun" $+ mhrun;
%
%ldraw  = lpath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "vd";
%odraws = opath $+ "\\" $+ lmodel $+ lprior $+ dataselstr $+ subTstr $+ "vds";
%
%open fhdraw = ^ldraw for read;
%
%nburn       = 101;        /* Number of initial draws to be discarded */
%hpdprob     = 0.90 ;
%
%drawrow = readr( fhdraw, nburn); 
%drawdim = cols(drawrow);
%drawrow = seekr( fhdraw, nburn);

drawmean = mean(hacdraw, 1);
drawstdd = std(hacdraw, 1);
drawcov = cov(hacdraw);

hpdprob = 0.90;
drawci = []

for i = 1:nvar
  drawci= [drawci, hpdint(hacdraw(:,i), hpdprob)];
end

%% Part 1: Compute the mean of x(i) and x^2(i)
%
%drawmean   = zeros(1,drawdim);
%drawsqmean = zeros(1,drawdim);
%drawcross  = zeros(drawdim,drawdim);
%
%ndraws = 0
%for i = 1:nblock
%
%   drawblock = ofhac[1:nuse+(i-1)*nuse,:];
%
%   drawmean  = drawmean    + sumc(drawblock)';
%   drawsqmean = drawsqmean + sumc(drawblock^2)';
%   drawcross  = drawcross  + drawblock'*drawblock;  
%
%   ndraws = ndraws + rows(drawblock);
%
%end
%
%drawmean   = drawmean/ndraws;
%drawsqmean = drawsqmean/ndraws;
%drawstdd   = sqrt(drawsqmean - (drawmean)^2);
%drawcov    = drawcross/ndraws - drawmean'*drawmean;




%/* Report Posterior Mean and stderror
%"MODEL       " lmodel;
%"PRIOR       " lprior;
%"MH-Run      " mhrun;
%"Subsample   " subT;
%
%"=============================================";
%" " ;
%
%load path = ^priopath para_names;
%"Parameter | Mean  |  StdD |";
%draw_names = { Y_c Y_i Y_q  Y_L Y_w Y_z Y_p Y_g Y_m 
%             Pi_c Pi_i Pi_q Pi_L Pi_w Pi_z Pi_P Pi_g Pi_m 
%               W_c W_i  W_q  W_L W_w W_z W_p W_g W_m 
%               I_c I_i  I_q  I_L I_w I_z I_p I_g I_m    
%               C_c C_i  C_q  C_L C_w C_z C_p C_g C_m 
%               R_c R_i  R_q  R_L R_w R_z R_p R_g R_m                       
%               L_c L_i  L_q  L_L L_w L_z L_p L_g L_m  };
%
%rows(drawmean);  cols(drawmean);
%
%outmat = draw_names'~ drawmean' ~ drawstdd' ;
%let mask[1,3] = 0 1 1;
%let fmt[3,3] =
%   "-*.*s " 8 4
%   "*.*lf " 8 4
%   "*.*lf " 8 4;
%
%d = printfm(outmat,mask,fmt);
%
%
%/* Part 2: HPD Interval
%*/
%cls;
%drawrow     = seekr( fhdraw, nburn); 
%drawci      = zeros(2,drawdim);
%
%j = 1;
%do until j > drawdim;
%
%   /* Read only the j'th column
%   */
%   drawrow   = seekr( fhdraw, nburn); 
%   drawblock = readr( fhdraw, 1);
%
%   drawcol = drawblock[1,j];
%   eofloop = 0;
%   ndraws  = 1;
%   
%   do until eofloop; 
%
%      drawblock = readr( fhdraw, nblock );
%      drawcol   = drawcol | drawblock[.,j];
%      ndraws    = ndraws + rows(drawblock);
%      eofloop   = eof(fhdraw);
%      locate 1,1;
%      "Part 2";
%      "ndraws" ndraws;
%      "Column" j;
%         
%   endo;
%   
%   drawci[.,j] = hpdint(drawcol,hpdprob);
%   j = j+1;
%
%endo;
%clear drawcol;
%
%/* Report Posterior Mean and stderror
%**
%*/
%cls;
%
%output file= c:\dge_sw\results\var_decomp_SW.txt reset;
%
%output on; 
%
%"MODEL       " lmodel;
%"PRIOR       " lprior;
%"MH-Run      " mhrun;
%"Subsample   " subT;
%"Coverage    " hpdprob;  
%"=============================================";
%" " ;
%
%load path = ^priopath para_names;
%"Parameter | Mean  |  StdD | CI(LOW) | CI(HIGH)";
%outmat = draw_names'~ drawmean' ~ real(drawstdd)' ~ drawci' ;
%let mask[1,5] = 0 1 1 1 1;
%let fmt[5,3] =
%   "-*.*s " 8 4
%   "*.*lf " 8 4
%   "*.*lf " 8 4
%   "*.*lf " 8 4
%   "*.*lf " 8 4;
%
%d = printfm(outmat,mask,fmt);
%
%create fhdraws = ^odraws with DRAWSUM, drawdim, 8;
%writer(fhdraws, real(drawmean | drawstdd | drawcov | drawci) );
%closeall fhdraws; 
%
%
%output off;
%
%end;
