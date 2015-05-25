function [numbersData] = IntronIntensities2(fillChannel,countChannel,cyclinChannel,intChannel);

%find amplitudes of transcription sites

contents = dir('data*');

m=1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        if obj.isGood == 0
            continue;
        end 
        
        if ~isfield(obj.metadata,'planeSpacing')
            disp('Need to tell me plane spacing. Please run function recordStackInfo.');
            return;
        end
        
        planeSpacing = obj.metadata.planeSpacing;
        
        %Find dapi mask
        dapiMask = obj.channels.dapi.processor.mask;
        
        %This is just the number of GAPDH mRNA for cell size
        numFillRNA = size(obj.channels.(fillChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
        numFillRNA = numFillRNA - inNuc;
        
        %Cyclin to determine cell cycle
        numCyclinRNA = size(obj.channels.(cyclinChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
        numCyclinRNA = numCyclinRNA - inNuc;
        
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
        numIntRNA = numel(find(a1&a2 == 1));
        
        if numIntRNA > 2 | numIntRNA == 0 %get rid of cells with >2 or 0 intron spots
            continue;
        end

        intronIntensities = obj.channels.(intChannel).metadata.gaussFitPostProc.amp(a1&a2);
        countIntensities = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(a0&a3);
        %countSig = obj.channels.(countChannel).metadata.gaussFitPostProc.sig(a3);
        avgIntron = mean(intronIntensities); %average intensity of intron probe
        
        avgRNA = mean(countIntensities); %average intensity of mRNA probe
        %txnSig = mean(countSig);
        
        numRnaPerTxnSite = countIntensities/singleSpotIntens;
        %num = numel(numRnaPerTxnSite)
        
        for k = 1:numel(numRnaPerTxnSite)
            %numbersData(m,:) = [i j obj.metadata.volumeRealUnits numFillRNA numRnaPerTxnSite(k)];
            numbersData(m,:) = [i j numFillRNA numCyclinRNA numRnaPerTxnSite(k)];
            m = m+1;
        end
        
    end
end

%figure; scatter(numbersData(:,3),numbersData(:,4));
%dlmwrite('GAPDH_Intens.txt',numbersData,'\t');