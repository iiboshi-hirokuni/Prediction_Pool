function Y = corrvc (M) 

[nrow, ncol] = size(M);

if nrow == ncol
  for i = 1:nrow
    for j = 1:ncol
      Y(i,j) = M(i,j)/sqrt(M(i,i))/sqrt(M(j,j));
    end
  end
else
  Y = M
end
