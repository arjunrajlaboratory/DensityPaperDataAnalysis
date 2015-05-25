function [] = extractalldata_heterokaryon(filename);

f = fopen(filename,'w');

% DO THIS IN TOP DIRECTORY BEFORE ANYTHING ELSE

% Cell type
spl = strsplit(pwd,'/');
cellType = char(spl(end));

mainDir = dir('1*');

fprintf(f,'exptType \t class \t date \t repNum \t dataNum \t objNum \t');
fprintf(f,'volume \t cytoGAPDH \t cytoGFP \t cytoGAS6 \n');

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
    
    spl = strsplit(pwd,'/');
    spl2 = strsplit(char(spl(end)),'_');
    class = char(spl2(end));
    date = char(spl2(1));
    
    if date == '131216'
        repNum = 1;
    else
        repNum = 2;
    end
    
    load(dataDir(1).name);
    
    fillChannel = 'nir';
    gfpChannel = 'tmr';
    gas6Channel = 'alexa';
    
    % DO THIS IN SUB DIRECTORY
    
    % DO THIS FOR EVERY CELL
    
    contents = dir('data*');
    cellNumber = 1;
    
    for k = 1:numel(contents)
        load(contents(k).name);
        for j = 1:numel(objects)
            
            obj = objects(j);
            
            if obj.isGood == 0
                continue;
            end
            
            %Find dapi mask
            dapiMask = obj.channels.dapi.processor.mask;
            if numel(dapiMask) == 0
                continue;
            end
            
            % Volume
            volume = obj.metadata.volumeRealUnits;
            
            % Total RNA
            totalRNA = obj.channels.(gfpChannel).numSpots;
            
            % Nuclear RNA
            nucRNA = numel(find(obj.channels.(gfpChannel).inNucleus(dapiMask)==1));
            
            % Cytoplasmic RNA
            cytoRNA = totalRNA - nucRNA;
            
            %GAS6
            totalGAS6 = obj.channels.(gas6Channel).numSpots;
            nucRNA = numel(find(obj.channels.(gas6Channel).inNucleus(dapiMask)==1));
            cytoGAS6 = totalGAS6 - nucRNA;
            
            % Total GAPDH
            totalGAPDH = obj.channels.(fillChannel).numSpots;
            
            % Nuclear GAPDH
            nucGAPDH = numel(find(obj.channels.(fillChannel).inNucleus(dapiMask)==1));
            
            % Cytoplasmic GAPDH
            cytoGAPDH = totalGAPDH - nucGAPDH;
            
            fprintf(f,'%s\t%s\t%s\t%d\t%d\t%d\t', cellType, class, date, repNum, k, j);
            fprintf(f,'%f\t%d\t%d\t%d\n', volume, cytoGAPDH, cytoRNA, cytoGAS6);
        end
        
    end
    cd('..');
end

fclose(f);