function [numbersData] = IntronFrequency2(fillChannel,countChannel,cyclinChannel,intChannel);

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
        
        mask = obj.object_mask.mask;
        
        %Find dapi mask
        dapiMask = obj.channels.dapi.processor.mask;
        
        %Only count spots outside the nucleus
        numFillRNA = size(obj.channels.(fillChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
        numFillRNA = numFillRNA - inNuc;
        
        %Intron spots are in the nucleus and colocalize with "count" spots
        %Returns indices of tmr spots in nucleus
        a1 = obj.channels.(intChannel).inNucleus(dapiMask);
        %Returns indices of tmr spots colocalized with alexa spots
        a2 = obj.channels.(intChannel).getColocalized(obj.channels.(countChannel),3,planeSpacing/.125);
        numIntRNA = numel(find(a1&a2 == 1));
        
        if numIntRNA > 2 %get rid of cells with >2 intron spots
            continue;
        end
                
        %numbersData(m,:) = [i j obj.metadata.volumeRealUnits numFillRNA numIntRNA];
        numbersData(m,:) = [i j numFillRNA numIntRNA];
        m = m+1;
    end
end

zeroSpots = numbersData(numbersData(:,4)==0,:);
oneSpot = numbersData(numbersData(:,4)==1,:);
twoSpots = numbersData(numbersData(:,4)==2,:);

mean0 = mean(zeroSpots(:,3));
mean1 = mean(oneSpot(:,3));
mean2 = mean(twoSpots(:,3));

er0 = 1.96*std(zeroSpots(:,3))/sqrt(numel(zeroSpots(:,3)));
er1 = 1.96*std(oneSpot(:,3))/sqrt(numel(oneSpot(:,3)));
er2 = std(twoSpots(:,3))/sqrt(numel(twoSpots(:,3)));

figure;errorbar([0 1 2], [mean0 mean1 mean2], [er0 er1 er2])