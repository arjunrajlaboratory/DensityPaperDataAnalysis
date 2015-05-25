function [traitCorr] = getTraitCorrTable(M,numIter,posSnpFreq,numGenesPerTrait)

traitMat = [];

for i = 1:numIter
    snpVec = rand(size(M,2),1)<posSnpFreq;
    expVec = M*snpVec;
    traitLevels = determineTraitLevels(expVec,numGenesPerTrait);
    traitMat = horzcat(traitMat,traitLevels');
end

traitCorr = corr(traitMat');