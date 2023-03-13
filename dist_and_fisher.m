% calculates correlation/partial correlation and fisher transforms the
% correlation coefficient
function [corrs]=dist_and_fisher(a,b,type,parcorrVec)

if nargin>3 % partial correltion
    corrs = partialcorr(a,b,parcorrVec);
else % regular correlation
    corrs=1-pdist2(a',b',type);
end

% fisher transformation
corrs=0.5*(log(1+corrs)-log(1-corrs));
end