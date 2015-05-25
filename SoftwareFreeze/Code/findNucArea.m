clear numbersData;
contents = dir('data*');
m=1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
       if obj.isGood == 0
            continue;
        end
        
        obj = objects(j);
        
        %Find dapi mask
        dapiMask = obj.channels.dapi.processor.mask;
        mask = obj.object_mask.mask;
        
        numbersData(m,:) = [obj.metadata.volumeRealUnits numel(find(mask)) numel(find(dapiMask))];
        m = m+1;
    end
end

dlmwrite('nuclear_area.txt',numbersData,'\t');