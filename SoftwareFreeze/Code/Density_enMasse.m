contents = dir('data*');

%fillChannel = 'nir';
%countChannel = 'alexa';
%intChannel = 'tmr';
%cyclinChannel = 'cy';

fillChannel = 'nir';
countChannel = 'tmr';
cyclinChannel = 'cy';


%Find the total number of objects in this folder
n = 0;
for i = 1:numel(contents)
    load(contents(i).name);
    n = n + numel(objects);
end

volumeData = zeros(n,6);

m = 1;

%Make a matrix with volume/number/density data
for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        %if (i==6 & j==2)
        %    continue;
        %end
        [vol numCount numFill dens] = findRNADensity(objects(j),fillChannel,countChannel);
        if objects(j).channels.(fillChannel).has_clearthreshold == 0
            numFill = 0;
        end
        if objects(j).channels.(countChannel).has_clearthreshold == 0
            numCount = 0;
        end
        %if objects(j).channels.(intChannel).has_clearthreshold == 0
        %    numInt = 0;
        %else
        %    numInt = numel(objects(j).channels.(intChannel).fitdataRNAonly.xp_fit);
        %end
        if objects(j).channels.(cyclinChannel).has_clearthreshold == 1
            numCyclin = numel(objects(j).channels.(cyclinChannel).fitdataRNAonly.xp_fit);
        else
            numCyclin = 0;
        end
        %volumeData(m,:)= [i j vol numFill numCount numInt numCyclin];
        volumeData(m,:)= [i j vol numFill numCount numCyclin];
        fprintf('Calculating density for object %d of %d\n',m,n);
        m = m+1;
    end
end

dlmwrite('densitymatrix_noBoundaryEnforce.txt',volumeData,'\t');