datafiles = dir('data*.mat');
color = 'cy';

m = 1;

clear intensMat

for i = 1:length(datafiles)
    fprintf('Processing file %d of %d\n',i,length(datafiles));
    load(datafiles(i).name);
    
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        
        img = obj.channelStk(color);
        
        %quantify intensity of each plane
        for k = 1:size(img,3)
            crop = img(:,:,k);
            nucCrop = crop.*uint16(obj.channels.cy.processor.mask);
            bkgrd = min(crop(find(crop(:))));
            crop = crop-bkgrd;
            
            dapiSum = sum(nucCrop(:));
            dapiAvg = dapiSum/numel(find(obj.channels.cy.processor.mask));
            
            vec(k) = dapiSum;
            %intensMat(m,k) = dapiSum;
        end
        vec = vec/max(vec);
        intensMat(m,:) = vec;
        m = m + 1;
    end
end

plot(intensMat');
xlabel('Image number','FontSize',20);
ylabel('Intensity','FontSize',20);
set(gca,'FontSize',20);