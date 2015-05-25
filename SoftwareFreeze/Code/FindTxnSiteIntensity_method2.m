% Locate z position for each transcription site, using intron channel

intronChannel = 'tmr';
cyclinChannel = 'cy';
fillChannel = 'nir';

contents = dir('data*');

m=1;
clear data

for i = 1:numel(contents) % each data file
    load(contents(i).name);
    for j = 1:numel(objects) % each image object
        obj = objects(j);
        
        if obj.isGood == 0
            continue;
        end
        
        dapiMask = obj.channels.dapi.processor.mask;
        
        if isfield(obj.metadata, 'volumeRealUnits')
            xVar = obj.metadata.volumeRealUnits;
        else
            xVar = obj.channels.(fillChannel).numSpots;
        end
        
        if ~isfield(obj.metadata,'transcriptionSitesManual')
            
            objects(j).metadata.transcriptionSitesManual.Xs = [];
            objects(j).metadata.transcriptionSitesManual.Ys = [];
            
            intensity = 0;
            numSites = 0;
            
            numCyclinRNA = size(obj.channels.(cyclinChannel).spotCoordinates,1);
            inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
            numCyclinRNA = numCyclinRNA - inNuc;
            
            data(m,:) = [i j intensity numSites numCyclinRNA xVar];
            m = m + 1;
            
            save(contents(i).name,'objects');
            continue;
        end
        im = obj.channelStk(intronChannel);
        Xs = obj.metadata.transcriptionSitesManual.Xs;
        Ys = obj.metadata.transcriptionSitesManual.Ys;
        xCoords = round(Xs);
        yCoords = round(Ys);
        
        numCyclinRNA = size(obj.channels.(cyclinChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
        numCyclinRNA = numCyclinRNA - inNuc;
        
        numSites = length(Xs);
        
        for k = 1:numel(obj.metadata.transcriptionSitesManual.Xs) % each transcription site
            x = xCoords(k);
            y = yCoords(k);
            
            imCropTiny = im((y-1):(y+1),(x-1):(x+1),:);
            maxList = [];
            for ii = 1:size(imCropTiny,3)
                maxList(ii) = max(max(imCropTiny(:,:,ii)));
            end
            [val plane] = max(maxList);
            
            imCropLarge = im((y-15):(y+15),(x-15):(x+15),plane);
            imCropLarge(6:25,6:25) = 0;
            background = median(imCropLarge(imCropLarge>0));
            
            intensity = val - background;
            
            data(m,:) = [i j intensity numSites numCyclinRNA xVar];
            m = m + 1;
        end
    end
end

dlmwrite('TxnSites_IntensityNumCyclin.txt',data,'\t');