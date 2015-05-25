function [numbersData] = getnumbers_novolume(fillChannel,cyclinChannel);

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
        % get GAPDH
        numFillRNA = size(obj.channels.(fillChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
        numFillRNA = numFillRNA - inNuc;
        
        % get Cyclin
        numCyclinRNA = size(obj.channels.(cyclinChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
        numCyclinRNA = numCyclinRNA - inNuc;
        
        % get transcription sites
        if isfield(obj.metadata,'transcriptionSitesManual')
            numTxnSites = length(obj.metadata.transcriptionSitesManual.Xs);
        else
            numTxnSites = 0;
        end
        
        numbersData(m,:) = [i j numFillRNA numCyclinRNA numTxnSites];
        m = m+1;
    end
end

dlmwrite('GapdhCyclinTxnSites.txt',numbersData,'\t');