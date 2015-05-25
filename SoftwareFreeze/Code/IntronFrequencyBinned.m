function [numbersData] = IntronFrequencyBinned(fillChannel,countChannel,cyclinChannel,intChannel,planeSpacing);

%clear numbersData;
contents = dir('data*');
%fillChannel = 'nir';
%countChannel = 'alexa';
%cyclinChannel = 'cy';
%intChannel = 'tmr';
%planeSpacing = .2;

m=1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
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
                
        numbersData(m,:) = [numFillRNA numIntRNA/2];
        m = m+1;
    end
end

[b,ix] = sort(numbersData(:,1));
sortx = numbersData(ix,:);
sz = size(sortx,1);
bound = ceil(sz/3);
small = sortx(1:bound,:);
large = sortx(end-bound+1:end,:);
avgSzSm = mean(small(:,1));
avgSzLg = mean(large(:,1));
avgFreqSm = mean(small(:,2));
avgFreqLg = mean(large(:,2));
stdFreqSm = std(small(:,2));
stdFreqLg = std(large(:,2));
strFreqSm = stdFreqSm/sqrt(size(small,1));
strFreqLg = stdFreqLg/sqrt(size(large,1));
figure; errorbar([avgSzSm avgSzLg],[avgFreqSm avgFreqLg],[strFreqSm strFreqLg]);
xlim([0 max(large(:,1))]);
ylim([0 1]);
xlabel('GAPDH','FontSize',20);
ylabel('Transcription probability per chromosome','FontSize',20);
%keyboard
%dlmwrite('IntronFrequencyBinned.txt',numbersData,'\t');