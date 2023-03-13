%/* filename:    CEEpm.g
%** description: The program maximizes the posterior density of the 
%**              DSGE model 
%** created:     10/05/2010
%*/

optflag = 2;   % 1:csminwel 2:Optimization Toolbox 3:Parallel


%******************************************************** 
% Starting Values for Maximization
%

iparaest = strcat(postpath,runname,'pm0.csv');
para = csvread(iparaest,0,1);

bounds = varargin.trspec(:,[2,3]);

ind1 = para > bounds(:,1);
ind2 = para < bounds(:,2);

ind1 = delif(ind1,varargin.pmask);
ind2 = delif(ind2,varargin.pmask);

if all(ind1) && all(ind2) == false
    disp('out of Bound');
    [  bounds(:,1) para  bounds(:,2) ]
    return;
end

%******************************************************** 
% Calculate posterior density at starting values

fprintf('Prior Density at Starting Value %f:\n',...
          priodens(para,pmean,pstdd,pshape));

[lnpY,retcode,obsmean,obsvar,shock] = evaldsge(para,YY,nshock,ZZ);

fprintf('Posterior at Starting Value: %f\n', lnpY);


figure(1000)  
  for j = 1:10
    subplot(4, 3, j)
        ti = 1981:0.25:1981+(nobs-1)/4;  
        plot(ti, obsmean(:,j),'k')       
  end
  
%*******************************************************
% Maximize the posterior density
%
tic;
if optflag == 1
    cc = 0;
    x0 = invtrans(para, trspec);  
    
    H0= eye(npara)*1E-1;
    nit = 1000;
    crit= 1E-4;
    
    [fh,xh,g,H,itct,fcount,csretcode] = csminwel('objfcn',x0,H0,[],crit,nit,varargin);

elseif optflag == 2
    x0 = invtrans(para, trspec);
    fun = @(x) objfcn(x,varargin);
    options = optimset('Display','iter','MaxIter', 10^5,'MaxFunEvals',10^6);
    [xh,fh,exitflag,output,g] = fminunc(fun,x0,options);
%    options = optimset('Display','iter','TolFun',1e-4, 'MaxIter', 10^5, 'MaxFunEvals', 10^6);
%    [xh,fh,exitflag,output] = fminsearch(fun,xh,options);
elseif optflag == 3
    x0 = invtrans(para, trspec);
    fun = @(x) objfcn(x,varargin);
    matlabpool open 2
    options = optimset('Display','iter','MaxIter', 10^5,'MaxFunEvals',10^6,'UseParallel','always');
    [xh,fh,exitflag,output,g] = fminunc(fun,x0,options);
    matlabpool close
end

%    options = optimset('Display','iter','MaxIter', 10^5,'MaxFunEvals',10^6,'DiffMaxChange',10.0);
%    [xh,fh,exitflag,output,g] = fminunc(fun,xh,options);

disp(sprintf('elapsed time is %f minutes',(toc/60)));

%*******************************************************
% Save parameter estimates
%

paraest = trans(xh, trspec);

oparaest = strcat(postpath,runname,'pm.csv');
csvwrite(oparaest, paraest);

minout(paraest,g,varargin,para_names,para);

delete H.dat
delete g1.mat
delete g2.mat
delete g3.mat

