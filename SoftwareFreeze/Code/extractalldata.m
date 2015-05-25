function [] = extractalldata(filename);

f = fopen(filename,'w');

% DO THIS IN TOP DIRECTORY BEFORE ANYTHING ELSE

% Cell type
spl = strsplit(pwd,'/');
cellType = char(spl(end));

mainDir = dir('1*');

fprintf(f,'cellType \t gene \t date \t dataNum \t objNum \t');
fprintf(f,'volume \t area \t totalGAPDH \t cytoGAPDH \t nucGAPDH \t');
fprintf(f,'nucArea \t nucIntensityTotal \t nucIntensityAvg \t');
fprintf(f,'numCyclin \t totalRNA \t cytoRNA \t nucRNA \t');
fprintf(f,'numTxnSites \t avgIntron \t avgExon \t numRnaPerTxnSiteAvgExon \n');

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
    spl = strsplit(pwd,'/');
    spl2 = strsplit(char(spl(end)),'_');
    gene = char(spl2(end));
    date = char(spl2(1));
    
    load(dataDir(1).name);
    
    if isfield(objects(1).channels,'alexa')
        if isfield(objects(1).channels,'nir')
            fillChannel = 'nir';
            countChannel = 'alexa';
            cyclinChannel = 'cy';
            intChannel = 'tmr';
        else
            fillChannel = 'alexa';
            countChannel = 'alexa';
            cyclinChannel = 'cy';
            intChannel = 'tmr';
        end
    else
        fillChannel = 'nir';
        countChannel = 'nir';
        cyclinChannel = 'cy';
        intChannel = 'tmr';
    end
    
    if numel(gene) == numel('POLR2A') & gene == 'POLR2A'
        fillChannel = 'nir';
        countChannel = 'tmr';
        cyclinChannel = 'cy';
        intChannel = [];
    end
    
    % DO THIS IN SUB DIRECTORY
    
    % DO THIS FOR EVERY CELL
    
    contents = dir('data*');
    cellNumber = 1;
    
    for k = 1:numel(contents)
        load(contents(k).name);
        for j = 1:numel(objects)
            
            obj = objects(j);
            
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
            
            % Volume
            if isfield(obj.metadata,'volumeRealUnits')
                volume = obj.metadata.volumeRealUnits;
            end
            
            % Cell area
            area = numel(find(mask))*.125*.125;
            
            % Cyclin count
            numCyclin = obj.channels.(cyclinChannel).numSpots;
            
            % Total RNA
            totalRNA = obj.channels.(countChannel).numSpots;
            
            % Nuclear RNA
            nucRNA = numel(find(obj.channels.(countChannel).inNucleus(dapiMask)==1));
            
            % Cytoplasmic RNA
            cytoRNA = totalRNA - nucRNA;
            
            % Total GAPDH
            totalGAPDH = obj.channels.(fillChannel).numSpots;
            
            % Nuclear GAPDH
            nucGAPDH = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
            
            % Cytoplasmic GAPDH
            cytoGAPDH = totalGAPDH - nucGAPDH;
            
            % Nuclear area
            nucArea = numel(find(dapiMask))*.125*.125;
            
            % Nuclear intensity
            [nucSum nucAvg] = getNucIntensity(obj,'dapi',fillChannel);
            
            if ischar(intChannel)
                
                %Number of transcription sites
                
                %Find the average intensity of a single mRNA
                notInNucIdx = find(obj.channels.(countChannel).inNucleus(dapiMask)==0);
                notInNucIntens = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(notInNucIdx);
                notInNucSig = obj.channels.(countChannel).metadata.gaussFitPostProc.sig(notInNucIdx);
                singleSpotSig = mean(notInNucSig);
                singleSpotIntens = mean(notInNucIntens);
                
                %Intron spots are in the nucleus and colocalize with "count" spots
                %Returns indices of tmr spots in nucleus
                a0 = obj.channels.(countChannel).inNucleus(dapiMask);
                a1 = obj.channels.(intChannel).inNucleus(dapiMask);
                %Returns indices of tmr spots colocalized with alexa spots
                a2 = obj.channels.(intChannel).getColocalized(obj.channels.(countChannel),3,planeSpacing/.125);
                a3 = obj.channels.(countChannel).getColocalized(obj.channels.(intChannel),3,planeSpacing/.125);
                
                numTxnSites = numel(find(a1&a2 == 1));
                
                intronIntensities = obj.channels.(intChannel).metadata.gaussFitPostProc.amp(a1&a2);
                countIntensities = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(a0&a3);
                avgIntron = mean(intronIntensities); %average intensity of intron probe
                if isnan(avgIntron)
                    avgIntron = 0;
                end
                
                avgRNA = mean(countIntensities); %average intensity of mRNA probe
                if isnan(avgRNA)
                    avgRNA = 0;
                end
                
                % Average number of RNA per transcription site
                numRnaPerTxnSiteAvgExon = avgRNA/singleSpotIntens;
                if isnan(numRnaPerTxnSiteAvgExon)
                    numRnaPerTxnSiteAvgExon = 0;
                end
                
            end
            
            fprintf(f,'%s\t%s\t%s\t%d\t%d\t', cellType, gene, date, k, j);
            if isfield(obj.metadata,'volumeRealUnits')
                fprintf(f,'%f\t%f\t%d\t%d\t%d\t', volume, area, totalGAPDH, cytoGAPDH, nucGAPDH);
            else
                fprintf(f,'%s\t%f\t%d\t%d\t%d\t', 'NA', area, totalGAPDH, cytoGAPDH, nucGAPDH);
            end
            fprintf(f,'%f\t%f\t%f\t', nucArea, nucSum, nucAvg);
            fprintf(f,'%d\t%d\t%d\t%d\t', numCyclin, totalRNA, cytoRNA, nucRNA);
            if ischar(intChannel)
                fprintf(f,'%d\t%f\t%f\t%f\n', numTxnSites, avgIntron, avgRNA, numRnaPerTxnSiteAvgExon);
            end
            
        end
    end
    
    cd('..');
end

fclose(f);