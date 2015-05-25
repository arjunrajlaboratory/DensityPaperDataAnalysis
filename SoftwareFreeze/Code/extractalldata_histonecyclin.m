function [] = extractalldata_histonecyclin(filename);

f = fopen(filename,'w');

% DO THIS IN TOP DIRECTORY BEFORE ANYTHING ELSE

% Cell type
spl = strsplit(pwd,'/');
cellType = char(spl(end));

mainDir = dir('1*');

fprintf(f,'exptType \t gene \t date \t fixative \t dataNum \t objNum \t');
fprintf(f,'numCyclin \t numHistone \n');

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
    % Fixative
    spl = strsplit(pwd,'/');
    spl2 = strsplit(char(spl(end)),'_');
    gene = char(spl2(3));
    date = char(spl2(1));
    fixative = char(spl2(2));
    
    cyclinChannel = 'cy';
    histoneChannel = 'alexa';
    
    
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
            
            % Cyclin count
            cyclin = obj.channels.(cyclinChannel).numSpots;
            inNuc = numel(find(obj.channels.(cyclinChannel).inNucleus(dapiMask)==1));
            numCyclin = cyclin - inNuc;
            
            % Histone count
            histone = obj.channels.(histoneChannel).numSpots;
            inNuc = numel(find(obj.channels.(histoneChannel).inNucleus(dapiMask)==1));
            numHistone = histone - inNuc;
            
            fprintf(f,'%s\t%s\t%s\t%s\t%d\t%d\t', cellType, gene, date, fixative, k, j);
            fprintf(f,'%d\t%d\n', numCyclin, numHistone);
        end
    end
    cd('..');
end

fclose(f);