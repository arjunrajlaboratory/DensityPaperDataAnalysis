% Get exon coordinates from intron coordinates

contents = dir('data*');

intChannel = 'tmr';
countChannel = 'alexa';
cyclinChannel = 'cy';
fillChannel = 'nir';

for i = 1:numel(contents)
    load(contents(i).name);
    disp(['Loading ' contents(i).name]);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        if obj.isGood == 0
            continue;
        end
        
        intIdx = obj.metadata.txnSiteIdx;
        
        ints = obj.channels.(intChannel).spotCoordinates(intIdx,1:2);
        
        exonsTotal = obj.channels.(countChannel).spotCoordinates(:,1:2);
        
        [nn,d] = knnsearch(exonsTotal,ints);
        
        nn = nn(d<3);
        
        objects(j).metadata.txnSiteIdx_Exons_Inferred = nn;
        
    end
    save(contents(i).name,'objects');
    disp(['Saving ' contents(i).name]);
end