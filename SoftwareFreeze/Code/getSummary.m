mainDir = dir('1*');

for i = 1:numel(mainDir) % For each large data directory
    
    if mainDir(i).isdir == 1
        disp(mainDir(i).name);
        cd(mainDir(i).name)
        
        numCol = 0;
        numDataFiles = 0;
        numObj = 0;
        numGood = 0;
        numProcessed = 0;
        numReviewed = 0;
        numVolume = 0;
        numRecInfo = 0;
        col = 0;
        
        filename = sprintf('%s_summary.txt',mainDir(i).name);
        f = fopen(filename,'w');
        
    else
        cd('..');
        continue;
    end
    
    dataDir = dir('data*');
    if numel(dataDir) < 1
        cd('..');
        continue;
    end
    %keyboard;
    for j = 1:numel(dataDir)
        load(dataDir(j).name)
        numDataFiles = numDataFiles + 1;
        if j == 1 & numel(objects) == 0
            start = 2;
            continue;
        elseif j == 1
            start = 1;
        end
        for k = 1:numel(objects)
            obj = objects(k);
            numObj = numObj + 1;
            if (j == start & k == 1)
                if isfield(obj.channels,'alexa')
                    fprintf(f,'alexa\t');
                    col = 'alexa';
                    numCol = numCol + 1;
                end
                if isfield(obj.channels,'cy')
                    fprintf(f,'cy\t');
                    if ~col
                        col = 'cy';
                    end
                    numCol = numCol + 1;
                end
                if isfield(obj.channels,'nir')
                    fprintf(f,'nir\t');
                    if ~col
                        col = 'nir';
                    end
                    numCol = numCol + 1;
                end
                if isfield(obj.channels,'cy')
                    fprintf(f,'tmr\t');
                    if ~col
                        col = 'tmr';
                    end
                    numCol = numCol + 1;
                end
                fprintf(f,'\n');
            end
            
            if obj.isGood == 1
                numGood = numGood + 1;
            end
            
            if obj.channels.(col).isProcessed == 1
                numProcessed = numProcessed + 1;
                if obj.channels.(col).processor.reviewed == 1
                    numReviewed = numReviewed + 1;
                end
            end
            
            if isfield(obj.metadata, 'volumeRealUnits')
                numVolume = numVolume + 1;
            end
            
            if isfield(obj.metadata, 'planeSpacing')
                numRecInfo = numRecInfo + 1;
            end

            
        end % end going through objects
        
    end % end going through data files
    
    fprintf(f, 'col = %s\n',col);
    fprintf(f, 'numDataFiles = %d\n', numDataFiles);
    fprintf(f, 'numObj = %d\n', numObj);
    fprintf(f, 'numGood = %d\n', numGood);
    fprintf(f, 'numProcessed = %d\n', numProcessed);
    fprintf(f, 'numReviewed = %d\n', numReviewed);
    fprintf(f, 'numVolume = %d\n', numVolume);
    fprintf(f, 'numRecInfo = %d\n', numRecInfo);
    
    fprintf('col = %s\n',col);
    fprintf('numCol = %d\n',numCol);
    fprintf('numDataFiles = %d\n', numDataFiles);
    fprintf('numObj = %d\n', numObj);
    fprintf('numGood = %d\n', numGood);
    fprintf('numProcessed = %d\n', numProcessed);
    fprintf('numReviewed = %d\n', numReviewed);
    fprintf('numVolume = %d\n', numVolume);
    fprintf('numRecInfo = %d\n', numRecInfo);
    fprintf('\n\n');
    
    fclose(f);
    cd('..');
end % end going through large data directories