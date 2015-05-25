function [fracBetweenParents] = getFractionOfTraitsThatFallBetweenParentsRnd(M,posSnpFreq,numGenesPerTrait)

parent1 = rand(size(M,2),1)<posSnpFreq;
parent2 = rand(size(M,2),1)<posSnpFreq;

parent1exp = M*parent1;
parent2exp = M*parent2;

parent1traits = determineTraitLevels(parent1exp,numGenesPerTrait);
parent2traits = determineTraitLevels(parent2exp,numGenesPerTrait);

for i = 1:length(parent1)
    if rand() > .5
        offspring(i) = parent1(i);
    else
        offspring(i) = parent2(i);
    end
end

offspring = offspring';


offspringexp = M*offspring;
offspringtraits = determineTraitLevels(offspringexp,numGenesPerTrait);

expCompare = horzcat(parent1exp,parent2exp,offspringexp,(parent1exp+parent2exp)/2);
traitCompare = horzcat(parent1traits',parent2traits',offspringtraits');

offspringBetweenParents = ((offspringtraits > parent1traits) & (offspringtraits < parent2traits)) | ...
    ((offspringtraits < parent1traits) & (offspringtraits > parent2traits));

fracBetweenParents = sum(offspringBetweenParents)/length(offspringBetweenParents);