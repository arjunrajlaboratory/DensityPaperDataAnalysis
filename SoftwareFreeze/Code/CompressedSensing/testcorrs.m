%% Matrix 1: 10 SNPs affect each gene, no overlap
% Determine matrix to use

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

% Multiply 1000 random vectors by this matrix, determine traits

% Low, medium, high frequencies of SNPs produce qualitatively similar
% correlations
lowFreq = getTraitCorrTable(M,1000,.1,4);
medFreq = getTraitCorrTable(M,1000,.5,4);
hiFreq = getTraitCorrTable(M,1000,.9,4);

%medFreq = getTraitCorrTableRandom(M,1000,.5,4);
v = randperm(length(medFreq));

%figure; imagesc(medFreq); colorbar;
figure; imagesc(medFreq(v,v)); colorbar;

% But the correlations are not identical
%figure; imagesc(hiFreq-medFreq); colorbar;

%% Matrix 2: 10 SNPs affect two genes. Groups of genes do not overlap.
% Again, qualitatively similar, but quantitatively different correlation
% profiles for low, medium, high frequency SNPs.
% There is significantly higher correlation between traits using this block
% system than above.

M2 = zeros(10,100);
M2(1:2,1:20) = rand(2,20);
M2(3:4,21:40) = rand(2,20);
M2(5:6,41:60) = rand(2,20);
M2(7:8,61:80) = rand(2,20);
M2(9:10,81:100) = rand(2,20);

lowFreq = getTraitCorrTable(M2,1000,.1,4);
medFreq = getTraitCorrTable(M2,1000,.5,3);
hiFreq = getTraitCorrTable(M2,1000,.9,4);

figure; imagesc(medFreq); colorbar;

%% Matrix 3: 10 SNPs affect 5 genes. Groups of genes do not overlap.
% Intense correlation.

M5 = zeros(10,100);
M5(1:5,1:50) = rand(5,50);
M5(6:10,51:100) = rand(5,50);

lowFreq = getTraitCorrTable(M5,1000,.1,4);
medFreq = getTraitCorrTable(M5,1000,.5,4);
hiFreq = getTraitCorrTable(M5,1000,.9,4);

figure; imagesc(lowFreq); colorbar;

%% Matrix 4: Still groups, but some overlap this time.
% Correlation pattern similar to matrix 3

Mol = zeros(10,100);
Mol(1:5,1:75) = rand(5,75);
Mol(6:10,26:100) = rand(5,75);

lowFreq = getTraitCorrTable(Mol,1000,.1,4);
medFreq = getTraitCorrTable(Mol,1000,.5,4);
hiFreq = getTraitCorrTable(Mol,1000,.9,4);

figure; imagesc(lowFreq); colorbar;

%% Matrix 5: Everything affects everything.
% Appears slightly more correlated than case (1), but less than the others.

Mmax = rand(10,100);

lowFreq = getTraitCorrTable(Mmax,1000,.1,4);
medFreq = getTraitCorrTable(Mmax,1000,.5,4);
hiFreq = getTraitCorrTable(Mmax,1000,.9,4);

figure; imagesc(medFreq); colorbar;