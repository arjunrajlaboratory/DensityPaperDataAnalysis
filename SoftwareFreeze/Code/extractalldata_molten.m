f = fopen('alldata.txt','w');

% DO THIS IN TOP DIRECTORY BEFORE ANYTHING ELSE

% Cell type
spl = strsplit(pwd,'/');
cellType = char(spl(end));

mainDir = dir('1*');

for i = 1:1%numel(mainDir)
    if mainDir(i).isdir == 1
        disp(mainDir(i).name);
        cd(mainDir(i).name)
        
        % Gene
        spl = strsplit(pwd,'/');
        spl2 = strsplit(char(spl(6)),'_');
        gene = char(spl2(2));
        
    else
        cd('..');
        continue;
    end
    dataDir = dir('data*');
    if numel(dataDir) < 1
        cd('..');
        continue;
    end
    
    load(dataDir(1).name);
    
    if isfield(objects(1).channels,'alexa')
        fillChannel = 'nir';
        countChannel = 'alexa';
        cyclinChannel = 'cy';
        intChannel = 'tmr';
    else
        fillChannel = 'nir';
        countChannel = 'nir';
        cyclinChannel = 'cy';
        intChannel = 'tmr';
    end
    
%     if ~isfield(objects(1).channels,'nir')
%         fillChannel = 'alexa';
%         if date == 130222
%             countChannel = 'tmr';
%             intChannel = [];
%         end
%     end
%     
%     if numel(gene) == numel('POLR2A') & gene == 'POLR2A'
%         fillChannel = 'nir';
%         countChannel = 'tmr';
%         cyclinChannel = 'cy';
%         intChannel = [];
%     end
    
    % DO THIS IN SUB DIRECTORY
    
    
    
    % Experiment date
    date = char(spl2(1));
    
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
                fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
                fprintf(f,'%s\t%f\n','volume',volume);
            end
            
            % Cell area
            area = numel(find(mask))*.125*.125;
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%f\n','area',area);
            
            % Cyclin count
            numCyclin = obj.channels.(cyclinChannel).numSpots;
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%d\n','numCyclin',numCyclin);
            
            % Cell cycle stage
            %if numCyclin > 100
            %    cellCycleStage = 'G2';
            %else
            %    cellCycleStage = 'G1';
            %end
            %fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            %fprintf(f,'%s\t%s\n','cellCycleStage',cellCycleStage);
            
            % Total RNA
            totalRNA = obj.channels.(countChannel).numSpots;
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%d\n','totalRNA',totalRNA);
            
            % Nuclear RNA
            nucRNA = numel(find(obj.channels.(countChannel).inNucleus(dapiMask)==1));
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%d\n','nucRNA',nucRNA);
            
            % Cytoplasmic RNA
            cytoRNA = totalRNA - nucRNA;
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%d\n','cytoRNA',cytoRNA);
            
            % Total GAPDH
            totalGAPDH = obj.channels.(fillChannel).numSpots;
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%d\n','totalGAPDH',totalGAPDH);
            
            % Nuclear GAPDH
            nucGAPDH = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%d\n','nucGAPDH',nucGAPDH);
            
            % Cytoplasmic GAPDH
            cytoGAPDH = totalGAPDH - nucGAPDH;
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%d\n','cytoGAPDH',cytoGAPDH);
            
            % Nuclear area
            nucArea = numel(find(dapiMask))*.125*.125;
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%f\n','nucArea',nucArea);
            
            % Nuclear intensity
            [nucSum nucAvg] = getNucIntensity(obj,'dapi',fillChannel);
            fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
            fprintf(f,'%s\t%f\n','nucSum',nucSum);
            fprintf(f,'%s\t%f\n','nucAvg',nucAvg);
            
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
                fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
                fprintf(f,'%s\t%d\n','numTxnSites',numTxnSites);
                
                intronIntensities = obj.channels.(intChannel).metadata.gaussFitPostProc.amp(a1&a2);
                countIntensities = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(a0&a3);
                avgIntron = mean(intronIntensities); %average intensity of intron probe
                
                avgRNA = mean(countIntensities); %average intensity of mRNA probe
                
                % Average intensity of transcription site (intron)
                fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
                fprintf(f,'%s\t%f\n','avgIntron',avgIntron);
                
                % Average intensity of transcription site (exon)
                fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
                fprintf(f,'%s\t%f\n','avgExon',avgRNA);
                
                % Average number of RNA per transcription site
                numRnaPerTxnSiteAvgExon = avgRNA/singleSpotIntens;
                fprintf(f,'%s\t%s\t%s\t%d\t',cellType,gene,date,cellNumber);
                fprintf(f,'%s\t%f\n','numRnaPerTxnSiteAvgExon',numRnaPerTxnSiteAvgExon);
                
            end
            
            cellNumber = cellNumber + 1;
        end
    end
    
    cd('..');
end

fclose(f);