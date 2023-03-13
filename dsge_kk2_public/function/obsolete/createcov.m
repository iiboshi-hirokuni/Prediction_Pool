function [omega] = createcov(para)

npara = max(size(para));
omega = zeros(npara, npara);
for i = 1:npara
  omega(i, i) = para(i)^2;
end
