posSnpFreq = .25;
numIter = 1000;
numGenesPerTrait = 4;
%M = normrnd(0,1,10,100);

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

% snpVec = rand(size(M,2),1)<posSnpFreq;
% 
 %normalExpLevels = [30; 40; 50; 60; 70; 80; 90; 100; 110; 120];
 normalExpLevels = [50; 100; 150; 200; 250; 300; 350; 400; 450; 500];
% 
% actualExpLevels = normalExpLevels + M*snpVec;

traitMat = [];

for i = 1:numIter
    snpVec = rand(size(M,2),1)<posSnpFreq;
    expVec = normalExpLevels + M*snpVec;
    traitLevels = determineTraitLevels(expVec,numGenesPerTrait);
    traitMat = horzcat(traitMat,traitLevels');
end

traitCorr = corr(traitMat');

figure; imagesc(traitCorr); colorbar;
