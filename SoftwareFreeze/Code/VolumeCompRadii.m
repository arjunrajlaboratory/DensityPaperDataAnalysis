function [numbersData] = VolumeCompEEF2(fillChannel,countChannel);

contents = dir('data*');

m=1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        if ~isfield(obj.metadata,'planeSpacing')
            disp('Need to tell me plane spacing. Please run function recordStackInfo.');
            return;
        end
        
        planeSpacing = obj.metadata.planeSpacing;
        
        if obj.isGood == 0
            continue;
        end
        
        mask = obj.object_mask.mask;
        
        %Find dapi mask
        dapiMask = obj.channels.dapi.processor.mask;
        
        %Only count spots outside the nucleus
        numFillRNA = size(obj.channels.(fillChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
        numFillRNA = numFillRNA - inNuc;
        
        numCountRNA = size(obj.channels.(countChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(countChannel).inNucleus(dapiMask)==1));
        numCountRNA = numCountRNA - inNuc;
        
        %vol = obj.metadata.volumeRealUnits;
        vol_r10 = obj.metadata.volumeRealUnits_r10;
        vol_r20 = obj.metadata.volumeRealUnits_r20;
        vol_r35 = obj.metadata.volumeRealUnits_r35;
        vol_r50 = obj.metadata.volumeRealUnits_r50;
                
        numbersData(m,:) = [vol_r10 vol_r20 vol_r35 vol_r50 numFillRNA];
        m = m+1;
    end
end

dlmwrite('VolumeComparison_Radii.txt',numbersData,'\t');