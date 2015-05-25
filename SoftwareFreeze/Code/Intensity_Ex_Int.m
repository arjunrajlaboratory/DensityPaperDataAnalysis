% Get exon coordinates from intron coordinates

contents = dir('data*');

intChannel = 'tmr';
countChannel = 'alexa';
cyclinChannel = 'cy';
fillChannel = 'nir';

m = 1; clear dat;

for i = 1:numel(contents)
    load(contents(i).name);
    disp(['Loading ' contents(i).name]);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        if obj.isGood == 0
            continue;
        end
        
        exIdx = obj.metadata.txnSiteIdx_Exons_Inferred;
        exs = obj.channels.(countChannel).spotCoordinates(exIdx,1:2);
        
        intsTotal = obj.channels.(intChannel).spotCoordinates(:,1:2);
        
        [nn,d] = knnsearch(intsTotal,exs);
        
        intensityEx = obj.channels.(countChannel).metadata.gaussFitPostProc.amp(exIdx);
        intensityInt = obj.channels.(intChannel).metadata.gaussFitPostProc.amp(nn);
        
        for k = 1:numel(nn)
            
            dat(m,:) = [i j k intensityEx(k) intensityInt(k) obj.channels.(fillChannel).numSpots obj.channels.(cyclinChannel).numSpots];
            m = m + 1;
        
        end
        
    end
end

dlmwrite('Intensity_Ex_Int.txt',dat,'\t');