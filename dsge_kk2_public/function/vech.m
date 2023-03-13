function Vector = vech(Matrix)
n = length(Matrix);
b = 0;
Vector = NaN(n*(n+1)/2,1);
for col = 1:n
    idx = transpose(1:col);
    Vector(b+idx) = Matrix((col-1)*n+idx);
    b = b+length(idx);
end
