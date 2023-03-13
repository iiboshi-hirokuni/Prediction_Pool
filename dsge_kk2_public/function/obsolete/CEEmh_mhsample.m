

iparaest = strcat(postpath,runname,'pm.csv');
paraest = csvread(iparaest);

imult = strcat(postpath,runname,'mu.csv');
sigscale = csvread(imult);

ihess = strcat(postpath,runname,'hs.csv');
hessian = csvread(ihess);

mkdir(strcat(resupath,runpath));

cc = 0.2; 
nsamples = nsim;

targetpdf = @(x) exp(mhstep_fcn2(x',varargin)); % –Œã•ª•z
proprnd = @(x) (x' + cc*(sigscale*randn(npara,1)))'; % ’ñˆÄ—”ƒTƒ“ƒvƒ‰

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

newnblock = nblock - lastGoodBlock;

tic;
for n = 1:newnblock
    [parasim, accept]= mhsample(paraest',nsamples,'pdf',targetpdf,'proprnd',proprnd,'symmetric',1);
    
    disp(sprintf('block %d of %d', (n+lastGoodBlock), nblock));
    disp(sprintf('Rejection rate:    %f\n',1.0-accept));
    fname = strcat(resupath, runname,'/-block--',num2str(n+lastGoodBlock));
    save(fname, 'parasim');
end

parasim_out = [];
likesim_out = [];
for n = (nburnin_block+1):nblock
    iname = strcat(resupath, runpath,'/-block--',num2str(n+lastGoodBlock));
    load(iname, 'parasim');
    parasim1 = parasim(mod(1:nsim, nthin)==0,:);
    parasim_out = [parasim_out; parasim1];
    likesim1 = [];
    for i = 1:size(parasim1,1)
        [ret1, ret2] = mhstep_fcn(parasim1(i,:)',varargin);
        likesim1(i,1) = ret2;
    end
    likesim_out = [likesim_out; likesim1];
end
disp(sprintf('elapsed time is %f minutes',(toc/60)));

oparasim = strcat(resupath, runpath, '/', runname, 'pa.csv');
util_csvwrite(oparasim, parasim_out, para_names);

olikesim = strcat(resupath, runpath, '/', runname, 'st.csv');
util_csvwrite(olikesim, likesim_out, {'postlike'});


%appendmode    = 0;
%lastGoodBlock = 0;
%
%if appendmode == 1
%    fname = strcat(runpath,'-block-',num2str(lastGoodBlock));
%    load(fname,'parasim');
%    paravals = parasim(nsim,1:npara);
%    % Set new starting value to last good parameter draw
%    for i = 1:npara
%        paraest(i) = paravals(i);
%    end   
%end
%
%% Run the posterior simulator
%
%mhstep(sigmult,hessian,nblock,nsim,cc0,cc,paraest,varargin,runname,...
%    appendmode,lastGoodBlock,resupath);
%
%parasim_out = [];
%likesim_out = [];
%for n = (nburnin_block+1):nblock
%    iname = strcat(resupath, runpath,'/-block-',num2str(n+lastGoodBlock));
%    load(iname, 'parasim', 'postsim', 'likesim', 'proppostsim', 'rej');
%    parasim1 = parasim(mod(1:nsim, nthin)==0,:);
%    likesim1 = likesim(mod(1:nsim, nthin)==0,:);
%    parasim_out = [parasim_out; parasim1];
%    likesim_out = [likesim_out; likesim1];
%end
%
%oparasim = strcat(resupath, runpath, '/', runname, 'pa.csv');
%util_csvwrite(oparasim, parasim_out, para_names);
%
%olikesim = strcat(resupath, runpath, '/', runname, 'st.csv');
%util_csvwrite(olikesim, likesim_out, {'postlike'});


