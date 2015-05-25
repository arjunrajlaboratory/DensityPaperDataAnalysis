%Summed dapi intensity

color1 = 'cy';
color2 = 'tmr';
m = 1;

contents = dir('data*');
dapiImgs = dir('dapi*.tif');
imgs1 = dir(sprintf('%s*.tif',color1));
imgs2 = dir(sprintf('%s*.tif',color2));

clear intensityData;

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        mask = obj.object_mask.mask;
        
        dapiImg = imread(dapiImgs(i).name);
        dapiImg = double(dapiImg);
        dapiImg = imcrop(dapiImg,obj.object_mask.boundingbox);
        dapiMask = maskWithDapi(dapiImg); 
        dapiMask(~mask)=0;
        
        area = numel(find(dapiMask==1));
        
        cellSize = numel(find(mask==1));
        
%         img1 = imread(imgs1(i).name);
%         img1 = double(img1);
%         img1 = imcrop(img1,obj.object_mask.boundingbox);
%         int1 = img1(dapiMask == 1);
%         bkg1 = img1(dapiMask == 0 & mask == 1);
%         int1sum = double(sum(int1(:)));
%         bkgrd1 = double(area*median(bkg1(:)));
%         total1 = int1sum-bkgrd1;
%         avg1 = total1/area;
        
        img2 = imread(imgs2(i).name);
        img2 = double(img2);
        img2 = imcrop(img2,obj.object_mask.boundingbox);
        int2 = img2(dapiMask == 1);
        bkg2 = img2(dapiMask == 0 & mask == 1);
        int2sum = sum(int2(:));
        bkgrd2 = area*median(bkg2(:));
        total2 = int2sum-bkgrd2;
        avg2 = total2/area;
        
        % 1)data num, 2)obj num, 3)area of cell, 4)summed total pol II intensity,
        % 5)summed phospho pol II intensity, 6) avg pol II, 7)avg phospho
        
        %intensityData(m,:) = [i j cellSize total1 total2 avg1 avg2];
        %intensityData(m,:) = [i j cellSize total1 avg1];
        intensityData(m,:) = [i j cellSize total2 avg2];
        m = m+1;
    end
end

dlmwrite('polIIintensity.txt',intensityData,'\t');

%figure; scatter(intensityData(:,3),intensityData(:,4),'filled');xlabel('cell area');ylabel('total pol II intensity');
%figure; scatter(intensityData(:,3),intensityData(:,5),'filled');xlabel('cell area');ylabel('total phospho intensity');
%pol2corr = corr2(intensityData(:,3),intensityData(:,4))
%phosphoCorr = corr2(intensityData(:,3),intensityData(:,5))