function [] = extractalldata_Click(filename);

f = fopen(filename,'w');

% DO THIS IN TOP DIRECTORY BEFORE ANYTHING ELSE

% Cell type
spl = strsplit(pwd,'/');
cellType = char(spl(end));

mainDir = dir('1*');

fprintf(f,'cellType \t date \t dataNum \t objNum \t');
fprintf(f,'volume \t area \t totalGAPDH \t cytoGAPDH \t nucGAPDH \t');
fprintf(f,'nucArea \t dapiIntensityTotal \t dapiIntensityAvg \t');
fprintf(f,'dtIntensityTotal \t dtIntensityAvg \t');
fprintf(f,'numCyclin \n');

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
    date = char(spl(end));
    
    
    load(dataDir(1).name);
    
    if numel(objects) == 0
        load(dataDir(2).name);
    end
    
    fillChannel = 'nir';
    clickChannel = 'tmr';
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
            else
                volume = 0;
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
            if totalGAPDH > 100
                [nucSum nucAvg] = getNucIntensity(obj,'dapi',fillChannel);
            else
                continue;
            end
            
            %Click intensity
            [clickSum clickAvg] = getNucIntensity_10percent(obj,clickChannel,fillChannel);
            
            
            fprintf(f,'%s\t%s\t%d\t%d\t', cellType, date, k, j);
            fprintf(f,'%f\t%f\t%d\t%d\t%d\t', volume, area, totalGAPDH, cytoGAPDH, nucGAPDH);
            fprintf(f,'%f\t%f\t%f\t', nucArea, nucSum, nucAvg);
            fprintf(f,'%f\t%f\t', clickSum, clickAvg);
            fprintf(f,'%d\n', numCyclin);
            
        end
    end
    
    cd('..');
end

fclose(f);