%function dapiopticalsliceNew
%for use with data from aTrous algorithm
%use minimum intensity pixel of entire image as background, and subtract that from intensity values
%for each object, identify 
%Within object mask, find the xy plane with the maximum sum of integrated
%DAPI intensity, and use the integrated DAPI intensity from that plane and store as column 3 of matrix DAPIintensities.

datafiles = dir('data*.mat');
color = 'alexa';

DAPIintensities=[];

globalobjcount=1;

for i = 1:length(datafiles)
    % Loop through all datafiles to process.
    fprintf('Processing file %d of %d\n',i,length(datafiles));
    load(datafiles(i).name);
    
    
    
    % Make sure there is at least one object.
    if ~isempty(objects)
        ob = objects(1);
        
        dapistk = readmm(ob.channels.(color).filename);
        
        dapistk.imagedata=dapistk.imagedata-min(dapistk.imagedata(:)); %subtract "background" (taken as minimum of entire image)
        
        for k = 1:numel(objects)
            
            obj=objects(k);
            
            % crop the stack to image_object bounding box
            objCrop = rectcropmulti(dapistk.imagedata,obj.object_mask.boundingbox); %crop
            
            dimensions=size(objCrop);
            
            % make values outside of segmentation border zero -- the mask is in 2D so need to first multiply by size of z dimension of bounding box
            %masked=objCrop.*uint16( repmat( objects(k).object_mask.mask, [1 1 dimensions(3)] ));
            masked=objCrop.*uint16( repmat( objects(k).channels.dapi.processor.mask, [1 1 dimensions(3)] ));
            
            % integrate intensity for each z-plane of the image
            xysum=sum(sum(masked,2));
            
            % find the plane with the greatest integrated intensity
            [maxDAPI,maxplane]=max(xysum);
            avgDAPI = maxDAPI/sum(objects(k).channels.dapi.processor.mask(:));
            maxplane

                      
            %fill out matrix for this data file
            DAPIintensities(globalobjcount,1)=globalobjcount; %record globalobjcount within data file
            DAPIintensities(globalobjcount,2)=i; %record datafile #
            DAPIintensities(globalobjcount,3)=maxDAPI;
            DAPIintensities(globalobjcount,4)=avgDAPI;
            %DAPIintensities(globalobjcount,5)=obj.metadata.volumebinary; %record volume
            %DAPIintensities(globalobjcount,6)=obj.channels.alexa.numSpots; %gapdh
                       
            
            %increase globalobjcount by 1 for next loop
            globalobjcount=globalobjcount+1;
            
         end;
        
    end;
    
end;

save('DAPIintensities','DAPIintensities');
