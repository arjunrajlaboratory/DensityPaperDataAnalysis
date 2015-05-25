function [numbersData] = getnumbers_novolume(fillChannel,countChannel,cyclinChannel,intChannel);

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
        if numel(dapiMask) == 0
            continue;
        end
        
        %Only count spots outside the nucleus
        numFillRNA = size(obj.channels.(fillChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
        numFillRNA = numFillRNA - inNuc;
        
        numCountRNA = size(obj.channels.(countChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(countChannel).inNucleus(dapiMask)==1));
        numCountRNA = numCountRNA - inNuc;
        
        numCyclinRNA = size(obj.channels.(cyclinChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
        numCyclinRNA = numCyclinRNA - inNuc;
        
        %Intron spots are in the nucleus and colocalize with "count" spots
        %Returns indices of tmr spots in nucleus
        a1 = obj.channels.(intChannel).inNucleus(dapiMask);
        %Returns indices of tmr spots colocalized with alexa spots
        a2 = obj.channels.(intChannel).getColocalized(obj.channels.(countChannel),3,planeSpacing/.125);
        numIntRNA = numel(find(a1&a2 == 1));
        
        numbersData(m,:) = [i j numFillRNA numCountRNA numCyclinRNA numIntRNA];
        m = m+1;
    end
end

dlmwrite('Numbers.txt',numbersData,'\t');