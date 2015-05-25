%% Everything affects everything, inherit in chunks

M = rand(10,100);
M = orth(M')';
posSnpFreq = .25;
numGenesPerTrait = 4;
numIter = 1000;

fracVec = zeros(1,numIter);

for i = 1:numIter
    fracVec(i) = getFractionOfTraitsThatFallBetweenParents(M,posSnpFreq,numGenesPerTrait);
end

plotTitle = sprintf('numGenes = %d, snpFreq = %f',numGenesPerTrait,posSnpFreq);

figure; hist(fracVec,50); title(plotTitle)

%% Everything affects everything, inherit randomly

M = rand(10,100);
M = orth(M')';
posSnpFreq = .25;
numGenesPerTrait = 4;
numIter = 1000;

fracVec = zeros(1,numIter);

for i = 1:numIter
    fracVec(i) = getFractionOfTraitsThatFallBetweenParentsRnd(M,posSnpFreq,numGenesPerTrait);
end

plotTitle = sprintf('numGenes = %d, snpFreq = %f',numGenesPerTrait,posSnpFreq);

figure; hist(fracVec,50); title(plotTitle)

%% A few SNPs affect each gene, no overlap; inherit in chunks

M = zeros(10,100);
M(1,1:10) = rand(1,10);
M(2,11:20) = rand(1,10);
M(3,21:30) = rand(1,10);
M(4,31:40) = rand(1,10);
M(5,41:50) = rand(1,10);
M(6,51:60) = rand(1,10);
M(7,61:70) = rand(1,10);
M(8,71:80) = rand(1,10);
M(9,81:90) = rand(1,10);
M(10,91:100) = rand(1,10);

M = orth(M')';

posSnpFreq = .25;
numGenesPerTrait = 4;
numIter = 1000;

fracVec = zeros(1,numIter);

for i = 1:numIter
    fracVec(i) = getFractionOfTraitsThatFallBetweenParents(M,posSnpFreq,numGenesPerTrait);
end

plotTitle = sprintf('numGenes = %d, snpFreq = %f',numGenesPerTrait,posSnpFreq);

figure; hist(fracVec,50); title(plotTitle)

%% A few SNPs affect each gene, no overlap; inherit randomly

M = zeros(10,100);
M(1,1:10) = rand(1,10);
M(2,11:20) = rand(1,10);
M(3,21:30) = rand(1,10);
M(4,31:40) = rand(1,10);
M(5,41:50) = rand(1,10);
M(6,51:60) = rand(1,10);
M(7,61:70) = rand(1,10);
M(8,71:80) = rand(1,10);
M(9,81:90) = rand(1,10);
M(10,91:100) = rand(1,10);

M = orth(M')';

posSnpFreq = .25;
numGenesPerTrait = 4;
numIter = 1000;

fracVec = zeros(1,numIter);

for i = 1:numIter
    fracVec(i) = getFractionOfTraitsThatFallBetweenParentsRnd(M,posSnpFreq,numGenesPerTrait);
end

plotTitle = sprintf('numGenes = %d, snpFreq = %f',numGenesPerTrait,posSnpFreq);

figure; hist(fracVec,50); title(plotTitle)