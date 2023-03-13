% creats vector of upper/lower triangular parts of a matrix
function [vec]=getTriangular(mat,Lower)
% mat: data matrix
% Lower: 1=lower triangolar, otherwise=upper triangualr

if Lower==1
    mat=tril(mat);
else
    mat=triu(mat);
end

mat(mat==0)=nan;
vec=mat(~isnan(mat));

end
