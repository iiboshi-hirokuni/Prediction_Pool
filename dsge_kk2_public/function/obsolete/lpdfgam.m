function ldens = lpdfgam(x,a,b)
% /* log GAMMA PDF
% */
% local ldens;
  ldens = -lngam(a) -a*log(b)+ (a-1)*log(x) -x/b ;
