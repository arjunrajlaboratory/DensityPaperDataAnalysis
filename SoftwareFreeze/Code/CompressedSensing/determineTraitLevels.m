function [traitLevels] = determineTraitLevels(vec,numGenesPerTrait)

% We have 10 genes
% If we assume the levels of 4 genes contribute to each trait, then we can
% have 210 independent traits. Let's start there for something concrete.

allCombos = nchoosek(1:length(vec),numGenesPerTrait);

traitLevels = zeros(1,size(allCombos,2));

for i = 1:size(allCombos,1)
    traitLevels(i) = sum(vec(allCombos(i,:)));
end

