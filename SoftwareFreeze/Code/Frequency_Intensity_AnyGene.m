function [freqPos,errorFreqPos,freqNeg,errorFreqNeg,intensPos,errorIntensPos,intensNeg,errorIntensNeg] = Frequency_Intensity_AnyGene(countChannel,cyclinChannel,intChannel);

% Intron Frequency

% Get cyclin, intron, and volume for each cell

contents = dir('data*');

m=1;
n=1;

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

        % Cyclin
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

        intronIntensities = obj.channels.(intChannel).metadata.gaussFitPostProc.amp(a1&a2);
        countIntensities = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(a0&a3);
        avgIntron = mean(intronIntensities); %average intensity of intron probe
        
        avgRNA = mean(countIntensities); %average intensity of mRNA probe
        
        numRnaPerTxnSite = countIntensities/singleSpotIntens;

        x(m,:) = [numCyclinRNA numIntRNA];
        m = m+1;
        
        if numIntRNA == 0 %get rid of cells with 0 intron spots for intensity
            continue;
        end
        
        for k = 1:numel(numRnaPerTxnSite)
            y(n,:) = [numCyclinRNA numRnaPerTxnSite(k)];
            n = n+1;
        end
    end
end

% Frequency
%keyboard;

cyclinPos = x(x(:,1)>100,:);
cyclinNeg = x(x(:,1)<100,:);

intPos = cyclinPos(:,2);
intPos(intPos>4)=[];
intNeg = cyclinNeg(:,2);
intNeg(intNeg>2)=[];

freqPos = mean(intPos);
errorFreqPos = std(intPos)/sqrt(numel(intPos));

freqNeg = mean(intNeg);
errorFreqNeg = std(intNeg)/sqrt(numel(intNeg));

% Intron Intensities

if ~exist('y')
    intensPos = 0;
    errorIntensPos = 0;
    intensNeg = 0;
    errorIntensNeg = 0;
else
    
    cyclinPos = y(y(:,1)>100,:);
    cyclinNeg = y(y(:,1)<100,:);
    
    intPos = cyclinPos(:,2);
    intPos(intPos>4)=[];
    intNeg = cyclinNeg(:,2);
    intNeg(intNeg>2)=[];
    
    intensPos = mean(intPos);
    errorIntensPos = std(intPos)/sqrt(numel(intPos));
    
    intensNeg = mean(intNeg);
    errorIntensNeg = std(intNeg)/sqrt(numel(intNeg));
end