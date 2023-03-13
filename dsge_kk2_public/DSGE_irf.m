%/*******************************************************/
%/*                                                     */
%/* Computing Impulse responses based on                */
%/* parameters of DSGE models                           */
%/*                                                     */
%/*******************************************************/

% description: The program converts the DSGE parameters
%              into impulse response functions

nirf = 20; 
nuse    = 100; % use every nuse observation

iparasim = strcat(resupath, runname, '/', runname, 'pa.csv');
fhpara = csvread(iparasim, 1, 0);

indseq0 = mod(1:(nblock-nburnin_block)*nsim, nuse);
indseq = (indseq0 ~= 0);
parasim_post  = delif(fhpara,indseq);
nasim = size(parasim_post,1);

hirf = zeros(nasim,nvar*nshock*nirf);
for i = 1:nasim
    % For each element of the block, compute the desired statistic
    para = parasim_post(i,1:npara)';
    ssirf = dsgeirf(para,nirf,nvar,nshock);
    hirf(i,:) = reshape(ssirf,nvar*nshock*nirf,1)';
end

irf_m0 = mean(hirf, 1);
irf_s0 = std(hirf, 1);

hpdprob_irf = 0.90;
drawci_irf = [];
size_irf = size(hirf,2);
for i = 1:size_irf
  drawci_irf= [drawci_irf, hpdint(hirf(:,i), hpdprob_irf)];
end

irf_hpd_l0 = drawci_irf(1,:);
irf_hpd_h0 = drawci_irf(2,:);

irf_m     = reshape(irf_m0,nirf,nvar*nshock);
irf_hpd_l = reshape(irf_hpd_l0,nirf,nvar*nshock);
irf_hpd_h = reshape(irf_hpd_h0,nirf,nvar*nshock);

%[irf_hpd_l(:,1),irf_m(:,1),irf_hpd_h(:,1)]

titlestr = {'Preference Shock: Z^b ' ,...
           'Exg-Goods Demand Shock: Z^g' ,...
		   'Wage Markup Shock: Z^w ' ,...
		   'Price Markup Shock: Z^p ' ,...
		   'Investment Price Shock: Z^{\nu}' ,...
		   'Monetary Policy Shock: Z^r' ,...
		   'Productivity Shock: Z^z ' ,... 
		   'Investment Shock: Z^{\psi}' ,...
           'Investment Price Shock: Z^{i}' ,...
		   'Financial Shock: Z^{efp}' ,...
		   'Net Worth Shock: Z^{nw}'};
%
ystr = {' -> Output' ,...
        ' -> Consumption' ,...
			' -> Investment' ,...  
			' -> Real Wage' ,...
			' -> Labor' ,...  
            ' -> Inflation' ,...
		    ' -> Investment Price' ,...
             ' -> Nominal Rate' ,...
             ' -> Real Borrowing', ...
             ' -> Loan Rate' };
             

for i = 1:nshock
   figure(3000+10*i)  
  for j = 1:nvar
    subplot(4, 3, j)
    plot(1:nirf, irf_m(:,nvar*(i-1)+j),'b',...
         1:nirf, irf_hpd_l(:,nvar*(i-1)+j),'-.r',...
         1:nirf, irf_hpd_h(:,nvar*(i-1)+j),'-.r')
    title(strcat(titlestr(i),ystr(j)))
  end
%   w = waitforbuttonpress;
end

