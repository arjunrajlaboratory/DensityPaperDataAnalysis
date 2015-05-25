contents = dir('data*');

intChannel = 'tmr';
countChannel = 'alexa';
cyclinChannel = 'cy';
fillChannel = 'nir';

%obj = objects(1);

m = 1;

for i = 1:numel(contents)
    load(contents(i).name);
    disp(['Loading ' contents(i).name]);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        if obj.isGood == 0
            continue;
        end
        
        if isfield(obj.metadata,'txnSiteIdx_Exon')
            continue;
        end
        
        im=imfuse(obj.channels.(intChannel).getImage,obj.channels.(countChannel).getImage);
        figure; imshow(im);
        hold on; scatter(obj.channels.(countChannel).spotCoordinates(:,2),obj.channels.(countChannel).spotCoordinates(:,1),'g');
        hold on; scatter(obj.channels.(intChannel).spotCoordinates(:,2),obj.channels.(intChannel).spotCoordinates(:,1),'r');
        
        zoom on;   % use mouse button to zoom in or out
        % Press Enter to get out of the zoom mode.
        % CurrentCharacter contains the most recent key which was pressed after opening
        % the figure, wait for the most recent key to become the return/enter key
        waitfor(gcf,'CurrentCharacter',13)
        
        [x,y] = ginput(4);
        close;
        % Locate the exon spots closest to x and y
        
        if numel(x) == 0
            exAvg = 0;
            exVal = [0 0 0 0];
            goodIdx = [];
        else
            
            exy = obj.channels.(countChannel).spotCoordinates(:,1);
            exx = obj.channels.(countChannel).spotCoordinates(:,2);
            
            exCoord = [exx exy];
            searchCoord = [x y];
            
            nn = knnsearch(exCoord,searchCoord);
            
            inNucIdx = find(obj.channels.(countChannel).inNucleus(obj.channels.dapi.processor.mask));
            
            overlap = ismember(nn,inNucIdx);
            goodIdx = nn(overlap);
            
            exV = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(goodIdx);
            
            exAvg = mean(exV);
            
            exVal = [0 0 0 0];
            exVal(1:numel(goodIdx)) = exV;
            
            exNum = numel(goodIdx);
            
            % Find intensity of regular mRNA spot
            notInNucIdx = find(obj.channels.(countChannel).inNucleus(obj.channels.dapi.processor.mask)==0);
            notInNucIntens = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(notInNucIdx);
            notInNucSig = obj.channels.(countChannel).metadata.gaussFitPostProc.sig(notInNucIdx);
            singleSpotSig = mean(notInNucSig);
            singleSpotIntens = mean(notInNucIntens);
            
        end
        
        objects(j).metadata.txnSiteIdx_Exon = goodIdx;
        
        %dat(m,:) = [i j intNum intAvg intVal obj.metadata.volumeRealUnits obj.channels.(cyclinChannel).numSpots];
        dat(m,:) = [i j exNum exAvg exVal singleSpotIntens obj.channels.(fillChannel).numSpots obj.channels.(cyclinChannel).numSpots];
        m = m + 1;
    end
    save(contents(i).name,'objects');
    disp(['Saving ' contents(i).name]);
end

dlmwrite('TxnSiteIntensityManual_Exon.txt',dat,'\t');