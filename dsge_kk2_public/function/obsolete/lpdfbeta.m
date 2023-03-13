function ldens = lpdfbeta(x,a,b)
% /* log BETA PDF
% */
% local ldens;
  ldens = lngam(a+b) - lngam(a) - lngam(b) + (a-1)*log(x) + (b-1)*log(1-x);
  