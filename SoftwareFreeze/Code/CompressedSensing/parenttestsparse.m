M = randn(250,1000);
M = 10*orth(M')';
posSnpFreq = .0385;

normalExpLevels = rand(250,1)*5+5;%%rand(250,1)*5;%ones(250,1);%[2; 3; 4; 5; 6; 7; 8; 9; 10; 11];
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

offspringBetweenParents = ((offspringexp > parent1exp) & (offspringexp < parent2exp)) | ...
    ((offspringexp < parent1exp) & (offspringexp > parent2exp));

%fracBetweenParents = sum(offspringBetweenParents)/length(offspringBetweenParents)

%% Now traits

Mtraits = randn(1000,250);
Mtraits = orth(Mtraits')';

parent1traits = Mtraits*parent1exp;
parent2traits = Mtraits*parent2exp;
offspringtraits = Mtraits*offspringexp;

corrWithParent1Traits = corr2(parent1traits,offspringtraits)
corrWithParent2Traits = corr2(parent2traits,offspringtraits)
corrWithMeanTraits = corr2((parent1traits+parent2traits)/2,offspringtraits)

%% Do these vectors even work for compressed sensing?

x = parent1geno;

figure; subplot(3,1,1); plot(x); title('original signal');

% measurement matrix
disp('Creating measurment matrix...');
A = randn(250,1000);
A = orth(A')';
disp('Done.');
	
% observations
y = A*x;

subplot(3,1,2); plot(y,'k'); title('measurements');


% initial guess = min energy
x0 = A'*y;

% solve the LP
tic
xp = l1eq_pd(x0, A, [], y, 1e-3);
toc


subplot(3,1,3); plot(xp,'r'); title('reconstructed signal');