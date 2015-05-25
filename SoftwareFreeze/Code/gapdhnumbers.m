contents = dir('data*');
m = 1;
clear gapdhData;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        gapdhData(m,:) = [obj.channels.nir.numSpots];
        m = m+1;
    end
end

%dlmwrite('gapdhnumbers.txt',gapdhData,'\t');

mu = mean(gapdhData);
sigma = std(gapdhData);

x = 1:max(gapdhData);
for i = 1:numel(x)
    y(i) = (1/sqrt(2*pi*sigma^2))*exp(-(x(i)-mu)^2/(2*sigma^2));
end

hist(gapdhData);hold on; scatter(x,y*30000);