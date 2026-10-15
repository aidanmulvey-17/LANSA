function [pval, meanR, R, medianR] = permTestMeanR(maps, nPerm)
%PERMTESTMEANR Conservation of laminar profiles across areas.
%   [pval, meanR, R, medianR] = permTestMeanR(maps, nPerm)
%   maps : depth x features x areas. For a single cell type, features = 1;
%          for the spectrolaminar motif, features = frequencies.
%   Each area's map is unrolled into a single vector and Pearson r is
%   computed between every pair of areas (R, areas x areas). Conservation is
%   summarized as the mean (and median) r over all pairs.
%
%   Null distribution: each area's map is circularly shifted along the depth
%   dimension by an independent random amount (keeping its internal
%   structure), and the mean pairwise r is recomputed nPerm times.
%   pval = (number of null values >= observed + 1) / (nPerm + 1).

[nDepth, nFeat, nAreas] = size(maps);
pairs = tril(true(nAreas), -1);

R = corr(reshape(maps, nDepth * nFeat, nAreas));
meanR   = mean(R(pairs));
medianR = median(R(pairs));

nullR = zeros(nPerm, 1);
for p = 1:nPerm
    shifted = maps;
    for a = 1:nAreas
        shifted(:, :, a) = circshift(maps(:, :, a), randi(nDepth) - 1, 1);
    end
    Rp = corr(reshape(shifted, nDepth * nFeat, nAreas));
    nullR(p) = mean(Rp(pairs));
end
pval = (sum(nullR >= meanR) + 1) / (nPerm + 1);
end
