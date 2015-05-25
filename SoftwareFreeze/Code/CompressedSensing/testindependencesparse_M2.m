M = zeros(250,1000);

j = 1;
for i = 1:250
    M(i,j:(j+3)) = randn(1,4);%ones(1,4);%randn(1,4);
    j = j+4;
end

M = orth(M')';

posSnpFreq = .0385;

normalExpLevels = rand(250,1)*5+1;

expMat = [];

for i = 1:1000
    snpVec = rand(size(M,2),1)<posSnpFreq;
    expVec = normalExpLevels + M*snpVec;
    expMat = horzcat(expMat,expVec);
end

traitCorr = corr(expMat');

figure; imagesc(traitCorr); colorbar;