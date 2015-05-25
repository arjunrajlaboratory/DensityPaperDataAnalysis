function [] = extractalldata_triptolide(filename)

f = fopen(filename,'w');

% DO THIS IN TOP DIRECTORY BEFORE ANYTHING ELSE

% Experiment type
spl = strsplit(pwd,'/');
drug = char(spl(end));

mainDir = dir('1*');

fprintf(f,'drug \t gene \t date \t conc \t dataNum \t objNum \t');
fprintf(f,'numTxnSites \t intensity \n');

for i = 1:numel(mainDir)
    if mainDir(i).isdir == 1
        disp(mainDir(i).name);
        cd(mainDir(i).name);
    else
        cd('..');
        continue;
    end
    
    spl = strsplit(pwd,'/');
    spl2 = strsplit(char(spl(end)),'_');
    date = char(spl2(1));
    gene = char(spl2(4));
    
    subDir = dir('1*');
    for sd = 1:numel(subDir)
        if subDir(sd).isdir == 1
            disp(subDir(sd).name);
            cd(subDir(sd).name);
        else
            %cd('..');
            continue;
        end
        
        spl = strsplit(pwd,'/');
        spl2 = strsplit(char(spl(end)),'_');
        conc = char(spl2(end));
        
        dataDir = dir('data*');
        if numel(dataDir) < 1
            cd('..');
            continue;
        end
        
        if ~strcmp(conc,'Ctrl') | strcmp(date,'141115')
            exonChannel = 'alexa';
            intronChannel = 'tmr';
            cyclinChannel = 'cy';
        else
            exonChannel = 'alexa';
            intronChannel = 'tmr';
            cyclinChannel = 'nir';
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
                if strcmp(date,'141105')
                    numCyclin = 0;
                else
                    cyclin = obj.channels.(cyclinChannel).numSpots;
                    inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
                    numCyclin = cyclin - inNuc;
                end
                
                im = obj.channelStk(intronChannel);
                Xs = obj.metadata.transcriptionSitesManual.Xs;
                Ys = obj.metadata.transcriptionSitesManual.Ys;
                
                if length(Xs) == 0
                    intensity = 0;
                    numTxnSites = 0;
                    
                    if strcmp(date,'141105')
                        numCyclin = '';
                    else
                        cyclin = obj.channels.(cyclinChannel).numSpots;
                        inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
                        numCyclin = cyclin - inNuc;
                    end
                    
                    fprintf(f,'%s\t%s\t%s\t%s\t%d\t%d\t', drug, gene, date, conc, k, j);
                    fprintf(f,'%d\t%f\n', numTxnSites, intensity);
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
                    
                    fprintf(f,'%s\t%s\t%s\t%s\t%d\t%d\t', drug, gene, date, conc, k, j);
                    fprintf(f,'%d\t%f\n', numTxnSites, intensity);
                end
            end
        end
        cd('..');
    end
    cd('..');
end

fclose(f);