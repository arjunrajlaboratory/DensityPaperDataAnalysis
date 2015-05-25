function [numbersData] = derek_getnumbers(fillChannel,countChannel,cyclinChannel,intChannel);

contents = dir('data*');

m=1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        if obj.isGood == 0
            continue;
        end
        
        mask = obj.object_mask.mask;
        
        %Find dapi mask
        dapiMask = obj.channels.dapi.processor.mask;
        if numel(dapiMask) == 0
            continue;
        end
        
        %Only count spots outside the nucleus
        numFillRNA = size(obj.channels.(fillChannel).spotCoordinates,1);
        
        numCountRNA = size(obj.channels.(countChannel).spotCoordinates,1);
        
        numCyclinRNA = size(obj.channels.(cyclinChannel).spotCoordinates,1);
        
        numIntRNA = size(obj.channels.(intChannel).spotCoordinates,1);
        
        numbersData(m,:) = [i j numFillRNA numCountRNA numCyclinRNA numIntRNA];
        m = m+1;
    end
end

dlmwrite('Numbers.txt',numbersData,'\t');