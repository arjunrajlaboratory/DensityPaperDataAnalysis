contents = dir('data*');

intChannel = 'tmr';
countChannel = 'alexa';
cyclinChannel = 'nir';
fillChannel = 'alexa';

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
        
        if isfield(obj.metadata,'txnSiteIdx')
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
        % Locate the intron spots closest to x and y
        
        if numel(x) == 0
            intNum = 0;
            intAvg = 0;
            intVal = [0 0 0 0];
            goodIdx = [];
        else
            
            inty = obj.channels.(intChannel).spotCoordinates(:,1);
            intx = obj.channels.(intChannel).spotCoordinates(:,2);
            
            intCoord = [intx inty];
            searchCoord = [x y];
            
            nn = knnsearch(intCoord,searchCoord);
            
            inNucIdx = find(obj.channels.(intChannel).inNucleus(obj.channels.dapi.processor.mask));
            
            overlap = ismember(nn,inNucIdx);
            goodIdx = nn(overlap);
            
            intV = obj.channels.tmr.metadata.gaussFitPostProc.amp(goodIdx);
            
            intAvg = mean(intV);
            
            intVal = [0 0 0 0];
            intVal(1:numel(goodIdx)) = intV;
            
            intNum = numel(goodIdx);
        end
            
            objects(j).metadata.txnSiteIdx = goodIdx;
            
            %dat(m,:) = [i j intNum intAvg intVal obj.metadata.volumeRealUnits obj.channels.(cyclinChannel).numSpots];
            dat(m,:) = [i j intNum intAvg intVal obj.channels.(fillChannel).numSpots obj.channels.(cyclinChannel).numSpots];
            m = m + 1;
        end
        save(contents(i).name,'objects');
        disp(['Saving ' contents(i).name]);
    end
    
    dlmwrite('TxnSiteIntensityManual.txt',dat,'\t');