%% Everything affects everything.

M = rand(10,100);
M = orth(M')';
posSnpFreq = .25;
numGenesPerTrait = 5;

parent1 = rand(size(M,2),1)<posSnpFreq;
parent2 = rand(size(M,2),1)<posSnpFreq;

parent1exp = M*parent1;
parent2exp = M*parent2;

parent1traits = determineTraitLevels(parent1exp,numGenesPerTrait);
parent2traits = determineTraitLevels(parent2exp,numGenesPerTrait);

%% Method 1: first half of SNPs from parent1, second half from parent2

offspring = vertcat(parent1(1:50),parent2(51:100));
offspringexp = M*offspring;
offspringtraits = determineTraitLevels(offspringexp,numGenesPerTrait);

expCompare = horzcat(parent1exp,parent2exp,offspringexp,(parent1exp+parent2exp)/2);
traitCompare = horzcat(parent1traits',parent2traits',offspringtraits');

traitCompare = traitCompare(randperm(size(traitCompare,1)),:);

dlmwrite('traitcompare_5050.txt',traitCompare,'delimiter','\t')

%% Method 2: randomly choose which SNPs go to offspring

for i = 1:length(parent1)
    if rand() > 0.5
        offspring(i) = parent1(i);
    else
        offspring(i) = parent2(i);
    end
end

offspringexp = M*offspring;
offspringtraits = determineTraitLevels(offspringexp,numGenesPerTrait);

expCompare = horzcat(parent1exp,parent2exp,offspringexp,(parent1exp+parent2exp)/2);
traitCompare = horzcat(parent1traits',parent2traits',offspringtraits');

dlmwrite('traitcompare_random.txt',traitCompare,'delimiter','\t')

%% New matrix. Only some SNPs affect each gene.

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

posSnpFreq = .25;

parent1 = rand(size(M,2),1)<posSnpFreq;
parent2 = rand(size(M,2),1)<posSnpFreq;

parent1exp = M*parent1;
parent2exp = M*parent2;

parent1traits = determineTraitLevels(parent1exp,numGenesPerTrait);
parent2traits = determineTraitLevels(parent2exp,numGenesPerTrait);

for i = 1:length(parent1)
    if rand() > 0.5
        offspring(i) = parent1(i);
    else
        offspring(i) = parent2(i);
    end
end

offspringexp = M*offspring;
offspringtraits = determineTraitLevels(offspringexp,numGenesPerTrait);

expCompare = horzcat(parent1exp,parent2exp,offspringexp,(parent1exp+parent2exp)/2);
traitCompare = horzcat(parent1traits',parent2traits',offspringtraits');

dlmwrite('traitcompare_m2_random.txt',traitCompare,'delimiter','\t')