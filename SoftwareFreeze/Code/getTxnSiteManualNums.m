contents = dir('data*');

intChannel = 'tmr';
countChannel = 'alexa';
cyclinChannel = 'nir';
fillChannel = 'alexa';

%obj = objects(1);

m = 1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        if obj.isGood == 0
            continue;
        end
        
        %Find dapi mask
        dapiMask = obj.channels.dapi.processor.mask;
        if numel(dapiMask) == 0
            continue;
        end
        
        goodIdx = obj.metadata.txnSiteIdx;
        
        if numel(goodIdx) == 0
            intAvg = 0;
            intVal = [0 0 0 0];
            intNum = 0;
        else
            intV = obj.channels.tmr.metadata.gaussFitPostProc.amp(goodIdx);
            
            intAvg = mean(intV);
            
            intVal = [0 0 0 0];
            intVal(1:numel(goodIdx)) = intV;
            
            intNum = numel(goodIdx);
        end
        
        numCountRNA = size(obj.channels.(countChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(countChannel).inNucleus(dapiMask)==1));
        numCountRNA = numCountRNA - inNuc;
        
        dat(m,:) = [i j intNum intAvg intVal obj.metadata.volumeRealUnits obj.channels.(cyclinChannel).numSpots];
        %dat(m,:) = [i j intNum intAvg intVal obj.channels.(fillChannel).numSpots obj.channels.(cyclinChannel).numSpots];
        %dat(m,:) = [i j intNum intAvg intVal obj.channels.(fillChannel).numSpots numCountRNA];
        m = m + 1;
    end
end

dlmwrite('TxnSiteIntensityManual.txt',dat,'\t');