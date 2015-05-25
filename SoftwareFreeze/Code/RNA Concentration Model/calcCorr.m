function [corr avgN] = calcCorr(expt)

% normalize
expt = expt ./ sum(sum(expt, 1));

% marginal probabilities of volume and mRNA
marV = sum(expt, 2)';
marN = sum(expt, 1);

% just the values of volume and mRNA (volume may not be the values input to
% the simulation as it is just the number of rows, but this does not matter
% for calculating the correlation coefficient)
valV = (1:size(expt, 1));
valN = 1:size(expt, 2);

% average values of volume and mRNA
avgV = marV * valV';
avgN = marN * valN';

% variances
varV = (valV - avgV).^2 * marV';
varN = (valN - avgN).^2 * marN';

% standard deviations
stdV = sqrt(varV);
stdN = sqrt(varN);

% this calculates the correlation between mRNA and volume
tmp = (valV - avgV)' * (valN - avgN);
corr = sum(sum(expt .* tmp)) / (stdV * stdN);

