contents = dir('data*');
fillChannel = 'nir';
countChannel = 'alexa';

%[datNum objNum volume numFillRNA numCountRNA]
volumeData = zeros(50,5);
n = 1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        fprintf('calculating volume for object %d\n',n);
        v = findVolUnbiased(objects(j),fillChannel);
        nFill = numel(objects(j).metadata.(fillChannel).x);
        nCount = numel(objects(j).metadata.(countChannel).x);
        volumeData(n,:) = [i j v nFill nCount];
        n = n + 1;
    end
end

dlmwrite('volume_unbiased_includeNucVol.txt',volumeData,'\t');