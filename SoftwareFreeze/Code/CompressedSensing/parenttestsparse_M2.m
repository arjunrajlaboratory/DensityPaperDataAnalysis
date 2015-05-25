M = zeros(250,1000);

j = 1;
for i = 1:250
    M(i,j:(j+3)) = randn(1,4);%ones(1,4);%randn(1,4);
    j = j+4;
end

M = orth(M')';

posSnpFreq = .0385;

normalExpLevels = rand(250,1)*5+1;%%rand(250,1)*5;%ones(250,1);%[2; 3; 4; 5; 6; 7; 8; 9; 10; 11];
%normalExpLevels = zeros(10,1);

parent1geno = rand(size(M,2),1)<posSnpFreq;
parent2geno = rand(size(M,2),1)<posSnpFreq;

parent1exp = normalExpLevels + M*parent1geno;
parent2exp = normalExpLevels + M*parent2geno;

offspringgeno = [];

%offspringgeno(1:500) = parent1geno(1:500);
%offspringgeno(501:1000) = parent2geno(501:1000);

for i = 1:length(parent1geno)
    if rand() > 0.5
        offspringgeno(i) = parent1geno(i);
    else
        offspringgeno(i) = parent2geno(i);
    end
end

offspringgeno = offspringgeno';

offspringexp = normalExpLevels + M*offspringgeno;

%figure; scatter((parent1exp+parent2exp)/2,offspringexp); refline(1,0);

%figure; scatter(parent1exp,offspringexp); refline(1,0);

corrWithParent1 = corr2(parent1exp,offspringexp)
corrWithParent2 = corr2(parent2exp,offspringexp)
corrWithMean = corr2((parent1exp+parent2exp)/2,offspringexp)
%hold on;
%scatter((parent1exp+parent2exp)/2,parent1exp,'r');
%scatter((parent1exp+parent2exp)/2,parent2exp,'m');

%figure; scatter(1:250,offspringexp);
%hold on;
%scatter(1:250,parent1exp,'r');
%scatter(1:250,parent2exp,'m');