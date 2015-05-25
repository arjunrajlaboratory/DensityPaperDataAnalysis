% 100 SNPs, 10 genes, 100 traits
% Traits need to be independent
% Linear: similar patterns of SNPs beget similar trait patterns

% Need to have some mapping between gene expression and phenotype/trait in
% order to be able to figure things out from these matrices... right?
% --- can still check linearity without that mapping

% Case 1: Redundancy. 10 SNPs change expression levels of each gene
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

Mtest = makeAndTestVector(M,size(M,2),3);

% Case 2: Groups of SNPs affect groups of genes, no overlap between groups
M2 = zeros(10,100);
M2(1:2,1:20) = rand(2,20);
M2(3:4,21:40) = rand(2,20);
M2(5:6,41:60) = rand(2,20);
M2(7:8,61:80) = rand(2,20);
M2(9:10,81:100) = rand(2,20);

M2test = makeAndTestVector(M2,size(M2,2),3);


M5 = zeros(10,100);
M5(1:5,1:50) = rand(5,50);
M5(6:10,51:100) = rand(5,50);

M5test = makeAndTestVector(M5,size(M5,2),3);


% Case 3: Each SNP affects groups of genes, groups can overlap (weak-ish 
% links between groups)
Mol = zeros(10,100);
Mol(1:5,1:75) = rand(5,75);
Mol(6:10,26:100) = rand(5,75);

Moltest = makeAndTestVector(Mol,size(Mol,2),3);

% Case 4: Everything affects everything
Mmax = rand(10,100);

Mmaxtest = makeAndTestVector(Mmax,size(Mmax,2),3);