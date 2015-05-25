function [] = extractalldata_freqintens_COL11A1_TUSC3(filename);

f = fopen(filename,'w');

% DO THIS IN TOP DIRECTORY BEFORE ANYTHING ELSE

% Cell type
spl = strsplit(pwd,'/');
cellType = char(spl(end));

mainDir = dir('1*');
mainDir = mainDir((end-2):end);

fprintf(f,'exptType \t gene \t date \t fixative \t dataNum \t objNum \t');
fprintf(f,'numTxnSites \t numCyclin \t xVar \t intensity \n');

for i = 1:numel(mainDir)
    if mainDir(i).isdir == 1
        disp(mainDir(i).name);
        cd(mainDir(i).name);
    else
        cd('..');
        continue;
    end
    dataDir = dir('data*');
    if numel(dataDir) < 1
        cd('..');
        continue;
    end
    
    % Experiment date
    % Gene
    % Fixative
    spl = strsplit(pwd,'/');
    spl2 = strsplit(char(spl(end)),'_');
    gene = char(spl2(3));
    date = char(spl2(1));
    fixative = char(spl2(2));
    
    if (length(gene) == length('EEF2') & gene == 'EEF2')
        exonChannel = 'alexa';
        intronChannel = 'tmr';
        cyclinChannel = 'nir';
        fillChannel = 'alexa';
    else
        exonChannel = 'alexa';
        intronChannel = 'tmr';
        cyclinChannel = 'cy';
        fillChannel = 'nir';
    end
    
    
    % DO THIS IN SUB DIRECTORY
    
    % DO THIS FOR EVERY CELL
    
    contents = dir('data*');
    cellNumber = 1;
    
    for k = 1:numel(contents)
        load(contents(k).name);
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
            
            % Cyclin count
            cyclin = obj.channels.(cyclinChannel).numSpots;
            inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
            numCyclin = cyclin - inNuc;
            
            % GAPDH or volume
            if ( (length(gene) == length('EEF2') & gene == 'EEF2') | (fixative == 'Form') )
                xVar = obj.metadata.volumeRealUnits;
            else
                totalGAPDH = obj.channels.(fillChannel).numSpots;
                nucGAPDH = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
                xVar = totalGAPDH-nucGAPDH;
            end
            
            im = obj.channelStk(intronChannel);
            Xs = obj.metadata.transcriptionSitesManual.Xs;
            Ys = obj.metadata.transcriptionSitesManual.Ys;
            
            if length(Xs) == 0
                intensity = 0;
                numTxnSites = 0;
                
                numCyclinRNA = size(obj.channels.(cyclinChannel).spotCoordinates,1);
                inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
                numCyclinRNA = numCyclinRNA - inNuc;
                
                fprintf(f,'%s\t%s\t%s\t%s\t%d\t%d\t', cellType, gene, date, fixative, k, j);
                fprintf(f,'%d\t%d\t%f\t%f\n', numTxnSites, numCyclin, xVar, intensity);
                continue;
            end
            
            xCoords = round(Xs);
            yCoords = round(Ys);
            
            numTxnSites = length(xCoords);
            
            % NOW FOR EACH TRANSCRIPTION SITE
            
            for jj = 1:length(xCoords)
                x = xCoords(jj);
                y = yCoords(jj);
                
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
                
                fprintf(f,'%s\t%s\t%s\t%s\t%d\t%d\t', cellType, gene, date, fixative, k, j);
                fprintf(f,'%d\t%d\t%f\t%f\n', numTxnSites, numCyclin, xVar, intensity);
            end
        end
    end
    cd('..');
end

fclose(f);