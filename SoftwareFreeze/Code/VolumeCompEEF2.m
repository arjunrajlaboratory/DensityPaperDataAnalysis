function [numbersData] = VolumeCompEEF2(fillChannel,countChannel);

contents = dir('data*');

m=1;

for i = 1:numel(contents)
    i
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
        
        vol = obj.metadata.volumeRealUnits;
        volEEF2 = obj.metadata.volumeRealUnits_alexa;
        volHalf = obj.metadata.volumeRealUnitsHalf;
        volExtra = obj.metadata.volumeRealUnitsExtraPts;
        volHalfExtra = obj.metadata.volumeRealUnitsHalfExtraPts;
                
        numbersData(m,:) = [vol volEEF2 volHalf volExtra volHalfExtra numFillRNA numCountRNA];
        m = m+1;
    end
end

dlmwrite('VolumeComparison.txt',numbersData,'\t');