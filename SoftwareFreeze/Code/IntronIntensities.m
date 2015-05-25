function [numbersData] = IntronIntensities(fillChannel,countChannel,cyclinChannel,intChannel);

%find amplitudes of transcription sites

contents = dir('data*');

m=1;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        if ~isfield(obj.metadata,'planeSpacing')
            disp('Need to tell me plane spacing. Please run function recordStackInfo.');
            return;
        end
        
        %Find dapi mask
        dapiMask = obj.channels.dapi.processor.mask;
        
        %This is just the number of GAPDH mRNA for cell size
        numFillRNA = size(obj.channels.(fillChannel).spotCoordinates,1);
        inNuc = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
        numFillRNA = numFillRNA - inNuc;
        
        %Find the average intensity of a single mRNA
        notInNucIdx = find(obj.channels.(countChannel).inNucleus(dapiMask)==0);
        notInNucIntens = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(notInNucIdx);
        notInNucSig = obj.channels.(countChannel).metadata.gaussFitPostProc.sig(notInNucIdx);
        singleSpotSig = mean(notInNucSig);
        singleSpotIntens = mean(notInNucIntens);
        
        %Intron spots are in the nucleus and colocalize with "count" spots
        %Returns indices of tmr spots in nucleus
        a1 = obj.channels.(intChannel).inNucleus(dapiMask);
        %Returns indices of tmr spots colocalized with alexa spots
        a2 = obj.channels.(intChannel).getColocalized(obj.channels.(countChannel),3,planeSpacing/.125);
        a3 = obj.channels.(countChannel).getColocalized(obj.channels.(intChannel),3,planeSpacing/.125);
        numIntRNA = numel(find(a1&a2 == 1));
        
        if numIntRNA > 2 | numIntRNA == 0 %get rid of cells with >2 or 0 intron spots
            continue;
        end
        
        intronIntensities = obj.channels.(intChannel).metadata.gaussFitPostProc.amp(a2);
        countIntensities = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(a3);
        countSig = obj.channels.(countChannel).metadata.gaussFitPostProc.sig(a3);
        avgIntron = mean(intronIntensities); %average intensity of intron probe
        
        avgRNA = mean(countIntensities); %average intensity of mRNA probe
        txnSig = mean(countSig);
        
        numRnaPerTxnSite = avgRNA/singleSpotIntens;
        
%         if numel(intronIntensities) == 0
%             avgIntron = 0;
%             avgRNA = 0;
%             numRnaPerTxnSite = 0;
%             txnSig = 0;
%         end
        
        numbersData(m,:) = [i j numFillRNA numRnaPerTxnSite avgIntron avgRNA singleSpotSig txnSig];
        m = m+1;
    end
end

x = numbersData(:,3:4);
[b,ix] = sort(x(:,1));
sortx = x(ix,:);
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
%keyboard
xlim([0 max(large(:,1))]);
ylim([0 max(large(:,2))]);
xlabel('GAPDH','FontSize',20);
ylabel('Txn site intensity (transcripts/site)','FontSize',20);

[a b c d e f] = strread(pwd, '%s %s %s %s %s %s', 'delimiter', '/');
[a1 b1 c1 d1 e1] = strread(char(f), '%s %s %s %s %s', 'delimiter', '_');

geneName = char(a(2,1));
date = char(a1);

%fileName = ['/Users/opadovan/Dropbox/Transcription Site Analysis/' geneName '_' date '_Intensity.png'];

%print(fileName,'-dpng');

%dlmwrite('IntronSpotIntensities.txt',numbersData,'\t');