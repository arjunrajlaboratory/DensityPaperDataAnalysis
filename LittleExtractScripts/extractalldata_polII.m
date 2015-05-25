function [] = extractalldata_polII(filename);

f = fopen(filename,'w');

% DO THIS IN TOP DIRECTORY BEFORE ANYTHING ELSE

% Cell type
spl = strsplit(pwd,'/');
spl2 = strsplit(char(spl(end)),'_');
cellType = char(spl2(2));
date = char(spl2(1));
gene = char(spl2(3));

fprintf(f,'cellType \t date \t dataNum \t objNum \t');
fprintf(f,'volume \t area \t totalGAPDH \t cytoGAPDH \t nucGAPDH \t');
fprintf(f,'nucArea \t dapiIntensityTotal \t dapiIntensityAvg \t');
fprintf(f,'totalRNA \t cytoRNA \t nucRNA \t');
fprintf(f,'numCyclin \n');


dataDir = dir('data*');

fillChannel = 'alexa';
countChannel = 'tmr';
cyclinChannel = 'cy';

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
        
        % Total RNA
        totalRNA = obj.channels.(countChannel).numSpots;
        
        % Nuclear RNA
        nucRNA = numel(find(obj.channels.(countChannel).inNucleus(dapiMask)==1));
        
        % Cytoplasmic RNA
        cytoRNA = totalRNA - nucRNA;
        
        
        fprintf(f,'%s\t%s\t%d\t%d\t', cellType, date, k, j);
        fprintf(f,'%f\t%f\t%d\t%d\t%d\t', volume, area, totalGAPDH, cytoGAPDH, nucGAPDH);
        fprintf(f,'%f\t%f\t%f\t', nucArea, nucSum, nucAvg);
        fprintf(f,'%f\t%f\t%f\t', totalRNA, cytoRNA, nucRNA);
        fprintf(f,'%d\n', numCyclin);
    end
end



fclose(f);