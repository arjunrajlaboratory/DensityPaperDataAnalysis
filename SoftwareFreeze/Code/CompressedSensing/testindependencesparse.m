M = randn(250,1000);
M = 10*orth(M')';
posSnpFreq = .0385;

normalExpLevels = rand(250,1)*5+5;

expMat = [];

for i = 1:1000
    snpVec = rand(size(M,2),1)<posSnpFreq;
    expVec = normalExpLevels + M*snpVec;
    expMat = horzcat(expMat,expVec);
end

traitCorr = corr(expMat');

figure; imagesc(traitCorr); colorbar;