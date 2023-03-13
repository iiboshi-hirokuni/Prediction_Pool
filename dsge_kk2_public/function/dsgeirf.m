function yyirfall = dsgeirf(para,nirf,nvar,nshock);
%/* INPUT:
%**   para   : parameter vector 
%**   nirf   : length of irf
%** OUTPUT:
%**   ssirf  : impulse response 
%**            the structural shocks.
%**            0 if retcode=0
%**            (nirf by 25)
%*/

global ZZ;

%/* Compute Covariance Matrix 
%** and its Chol decomposition
%*/ 
sig_chol = eye(nshock); 
sig_chol(nshock,nshock) = 1/4;%  monetary

% solve DSGE

[T1,TC,T0,RC] = dsgesolv(para);

% ZZ = zeros(nvar, size(T1, 1));
% ZZ(1,1) = 1;               % output gap
% ZZ(2,2) = 4;               % inflation
% ZZ(3,3) = 1;               % wage
% ZZ(4,6) = 1/10;            % investment
% ZZ(5,7) = 1;               % consumption
% ZZ(6,8) = 4;               % nominal rate
% ZZ(7,10) = 1;              % labor


yyirfall = zeros(nirf,nvar*nshock); 

for sh_ind = 1:nshock
    impact = sig_chol(:,sh_ind);
    yyirf  = zeros(nirf,nvar);
    dyyirf  = zeros(nirf,nvar);
    ss = T0*impact;
    s  = [ss;zeros(size(ss,1),1)];
    dyyirf(1,:) = ZZ*s;
    yyirf(1,:) = dyyirf(1,:);
    
    for t = 2:nirf
        ss1 = T1*ss;     
        s = [ss1;ss];
        dyyirf(t,:) = (ZZ*s)';
        yyirf(t,1:4) = yyirf(t-1,1:4)+dyyirf(t,1:4);
        yyirf(t,5:10) = dyyirf(t,5:10);
        ss = ss1;
    end
    
    yyirfall(:,1+nvar*(sh_ind-1):nvar*sh_ind) = yyirf;
    
%     for j = 1:9
%       subplot(3, 3, j)
%       plot(1:nirf, yyirf(:,j),'b')
%       title(strcat(titlestr(sh_ind),ystr(j)))
%     end
%   w = waitforbuttonpress;
end

% 
% for sh_ind = 1:nshock
%    
%     impact = sig_chol(:,sh_ind);
%     yyirf  = zeros(nirf,nvar);
%     ss = T0*impact;
%     yyirf(1,:) = ZZ*ss;
%     
%     for t = 2:nirf
%         ss = T1*ss;     
%         yyirf(t,:) = (ZZ*ss)';
%     end
%     yyirfall(:,1+nvar*(sh_ind-1):nvar*sh_ind) = yyirf;
% end

