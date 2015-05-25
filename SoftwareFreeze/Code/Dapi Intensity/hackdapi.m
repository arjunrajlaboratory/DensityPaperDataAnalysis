%function dapiopticalsliceNew
%for use with data from aTrous algorithm
%use minimum intensity pixel of entire image as background, and subtract that from intensity values
%for each object, identify 
%Within object mask, find the xy plane with the maximum sum of integrated
%DAPI intensity, and use the integrated DAPI intensity from that plane and store as column 3 of matrix DAPIintensities.

datafiles = dir('data*.mat');
color = 'alexa';
gapdhColor = 'nir';

%DAPIintensities=[];
clear DAPIintensities

spl = strsplit(pwd,'/');
date = str2num(char(spl(end)));

globalobjcount=1;

for i = 1:length(datafiles)
    % Loop through all datafiles to process.
    fprintf('Processing file %d of %d\n',i,length(datafiles));
    load(datafiles(i).name);
    
    
    
    % Make sure there is at least one object.
%     if ~isempty(objects)
%         ob = objects(1);
%         
%         dapistk = readmm(ob.channels.(color).filename);
%         
%         dapistk.imagedata=dapistk.imagedata-min(dapistk.imagedata(:)); %subtract "background" (taken as minimum of entire image)
        
        for k = 1:numel(objects)
            
            obj=objects(k);
            
            if ~obj.isGood
                continue;
            end
            
            minPlane = min(obj.channels.(gapdhColor).spotCoordinates(:,3));
            maxPlane = max(obj.channels.(gapdhColor).spotCoordinates(:,3));
            cropPlane = round((maxPlane-minPlane)/2)
            %cropPlane = 6;
            
            if cropPlane == 0
                continue;
            end
            
            img=obj.channelStk(color);
            crop=img(:,:,cropPlane); 
            crop=crop.*uint16(obj.object_mask.mask);
            
            bkgrd = min(crop(find(crop(:))));
            
            crop = crop-bkgrd;
            
            nucCrop = crop.*uint16(obj.channels.dapi.processor.mask);
            
            dapiSum = sum(nucCrop(:));
            dapiAvg = dapiSum/numel(find(obj.channels.dapi.processor.mask));
      
            %fill out matrix for this data file
            %DAPIintensities(globalobjcount,1)=globalobjcount; %record globalobjcount within data file
            %DAPIintensities(globalobjcount,2)=i; %record datafile #
            %DAPIintensities(globalobjcount,1)=obj.metadata.volumeRealUnits;
            %DAPIintensities(globalobjcount,1)=date;
            DAPIintensities(globalobjcount,1)=date;
            DAPIintensities(globalobjcount,2)=obj.channels.nir.numSpots;
            DAPIintensities(globalobjcount,3)=dapiSum;
            DAPIintensities(globalobjcount,4)=dapiAvg;
            
            
            %increase globalobjcount by 1 for next loop
            globalobjcount=globalobjcount+1;
            
         end;
        
    %end;
    
end;
%DAPIintensities(:,1)=date;
name = [upper(color) 'intensities'];
name2 = ['sum_avg_gapdh_' color '.txt'];
save(name,'DAPIintensities');
dlmwrite(name2,DAPIintensities,'precision','%.6f')
