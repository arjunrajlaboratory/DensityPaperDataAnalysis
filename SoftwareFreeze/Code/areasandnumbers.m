clear;
contents = dir('data*');
%fillChannel = 'nir';
countChannel = 'alexa';
%countChannel = 'tmr';
%cyclinChannel = 'cy';
intChannel = 'tmr';

m=1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        mask = obj.object_mask.mask;
        
        area = numel(find(mask));
        
        numbersData(m,:) = [area obj.channels.(countChannel).numSpots obj.channels.(intChannel).numSpots];
        m = m+1;
    end
end

dlmwrite('Area_Alexa_Tmr.txt',numbersData,'\t');
%dlmwrite('areamatrix_excludeNucleus.txt',areaData,'\t');