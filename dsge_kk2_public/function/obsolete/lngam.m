function lng = lngam(x)
% /* This procedure computes the log-gamma function
% */
% local k,lng;
%   if x <= 3;
%      lng= log(gamma(x));
%   else
%      k = floor(x-2);
%      lng= sum(log(x-1:-1:k)) + log(gamma(x-k));
%   end
lng = gammaln(x);
