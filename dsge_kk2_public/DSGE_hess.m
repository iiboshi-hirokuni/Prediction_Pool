%----------------------------------------------------------------------
% description: The program computes the hessian at the posterior mode
%----------------------------------------------------------------------

iparaest = strcat(postpath,runname,'pm.csv');
paraest = csvread(iparaest);

parflag = 0;

%---- mhcov.m ----
%function [sigmult,hmat,penalt] = mhcov(para,vs)

tic;
if parflag == 0
    [hessian,errorij] = hessn('hessn_fcn',paraest,varargin);
elseif parflag == 1
    [hessian,errorij] = hessn_para('hessn_fcn',paraest,varargin);
end
disp(sprintf('elapsed time is %f minutes',(toc/60)));

X = -hessian;

[u,s,v] = svd(X,0);

rankHHM = sum(varargin.pmaskinv);

invHHMdet = 1;

for i = 1:length(paraest)
    if i > rankHHM
        s(i,i) = 0;
    else
        s(i,i) = 1.0/s(i,i);
        invHHMdet = invHHMdet*s(i,i);
    end
end

invHHM = u*s*u';
sigmult = u*sqrt(s);
penalt = (rankHHM/2)*log(2*pi) + 0.5*log(invHHMdet);

disp('Determinant of variance matrix');
disp(invHHMdet);
disp('s.e. at posterior mode');
disp(sqrt(diag(invHHM)));
disp('Post Mode Penalty');
disp(penalt);

%-------------------

omult = strcat(postpath,runname,'mu.csv');
csvwrite(omult, sigmult);

ohess = strcat(postpath,runname,'hs.csv');
csvwrite(ohess, hessian);

%save(hessname,'sigmult','hessian','penalt');

