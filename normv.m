function nv = normv(v)
% v is m x n matrix. Each row represent a vector to be normalized
nv = v ./ repmat(sqrt(sum(v.^2, 2)), 1, size(v, 2));
end