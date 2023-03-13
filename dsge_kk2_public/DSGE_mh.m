disp('=============================================================');
disp('   Start MH algorithm of DSGE with Frictions');
disp('=============================================================');



if use_post_mode == 1
   iparaest = strcat(postpath,runname,'pm.csv');
   paraest = csvread(iparaest);
else
   iparaest = strcat(postpath,runname,'pm0.csv');
   paraest = csvread(iparaest);
end
   
if use_post_hess == 1
   imult = strcat(postpath,runname,'mu.csv');
   sigmult = csvread(imult);
else
   sigmult = 0.05*diag(abs(paraest));
   
end


ihess = strcat(postpath,runname,'hs.csv');
hessian = csvread(ihess);

mkdir(strcat(resupath,runpath));

%load(hessname,'sigmult','hessian','penalt');

cc0    = cc; %0.25;  %% adjustment coefficient of Random Walk MH algorithm

% Code block to restart the MHRW if it is interrupted or wanders
% off into a bad region. If a resume, set appendmode to 1 and set the
% value for the last good block

appendmode    = 0;
lastGoodBlock = 0;

if appendmode == 1
    fname = strcat(runpath,'-block-',num2str(lastGoodBlock));
    load(fname,'parasim');
    paravals = parasim(nsim,1:npara);
    % Set new starting value to last good parameter draw
    for i = 1:npara
        paraest(i) = paravals(i);
    end   
end

% Run the posterior simulator

tic;
mhstep(sigmult,hessian,nblock,nburnin_block,nsim,cc0,cc,paraest,varargin,runpath,...
    appendmode,lastGoodBlock,resupath);

parasim_out = [];
likesim_out = [];
% leva_out = []; 

for n = (nburnin_block+1):nblock
    iname = strcat(resupath, runpath,'/-block-',num2str(n+lastGoodBlock));
    load(iname, 'parasim', 'postsim', 'likesim', 'proppostsim', 'rej');
    parasim1 = parasim(mod(1:nsim, nthin)==0,:);
    postsim1 = postsim(mod(1:nsim, nthin)==0,:);
    likesim1 = likesim(mod(1:nsim, nthin)==0,:);
%     levasim1 = levasim(mod(1:nsim, nthin)==0,:);
    parasim_out = [parasim_out; parasim1];
    likesim_out = [likesim_out; postsim1 likesim1];
end

oparasim = strcat(resupath, runpath, '/', runname, 'pa.csv');
util_csvwrite(oparasim, parasim_out, para_names);

est_date = datestr(date);   
oparasim1 = strcat('c:\DSGE\DSGE_KK/results/', runname,'_',num2str(nsim*(nblock-nburnin_block)),'_',est_date,'pa.csv');
util_csvwrite(oparasim1, parasim_out, para_names);

olikesim = strcat(resupath, runpath, '/', runname, 'st.csv');
util_csvwrite(olikesim, likesim_out, {'postlike'});
disp(sprintf('elapsed time is %f minutes',(toc/60)));

