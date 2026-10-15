function X_full = buildFeatureMatrix(raw_anatomy, ref)
%BUILDFEATUREMATRIX PLSR predictors from laminar cell densities.
%   X_full = buildFeatureMatrix(raw_anatomy) returns 11 predictors for a
%   [bins x 4] density matrix with columns PV, CB, CR, NG (pyramidal):
%      1-4   z-scored PV, CB, CR and NG densities
%      5-10  interactions: products of the z-scored densities
%            (PV×CB, PV×CR, PV×NG, CB×CR, CB×NG, CR×NG)
%      11    z-scored E/I balance, NG / (PV + CB + CR)
%   Z-scores are computed across laminar bins.
%
%   X_full = buildFeatureMatrix(raw_anatomy, ref) z-scores with fixed
%   reference statistics instead (struct with fields mu, sigma [1 x 4],
%   mu_ei, sigma_ei), so that control and perturbed densities are on the
%   same scale (schizophrenia analysis).

inhibs = sum(raw_anatomy(:, 1:3), 2);
pyrams = raw_anatomy(:, 4);
ei = (pyrams + eps) ./ inhibs;

if nargin < 2
    X_z = zscore(raw_anatomy);
    EI  = zscore(ei);
else
    X_z = (raw_anatomy - ref.mu) ./ ref.sigma;
    EI  = (ei - ref.mu_ei) ./ ref.sigma_ei;
end

PV = X_z(:, 1);  CB = X_z(:, 2);  CR = X_z(:, 3);  NG = X_z(:, 4);
ints = [PV.*CB, PV.*CR, PV.*NG, CB.*CR, CB.*NG, CR.*NG];

X_full = [X_z, ints, EI];
end
