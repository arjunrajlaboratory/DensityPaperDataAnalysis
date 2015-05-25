color = 'tmr';
contents = dir('data*');
m = 1;
clear intData;


for i = 1:numel(contents)
    load(contents(i).name);
    fprintf('loading %s\n',contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        
        img = obj.channels.(color).maxlaplaceimage;
        imagesc(img);
        hold on;
        
        %Find dapi mask
        dapiStk = obj.channelStk('dapi');
        dapiStk = max(dapiStk,[],3);
        dapiMask1 = maskWithDapi(dapiStk);
        dapiMask = bwareaopen(dapiMask1,2000);
        dapiMask = imresize(dapiMask,expFactor);
        
        bw = bwperim(dapiMask);
        [x,y] = find(bw==1);
        scatter(y,x,'r');
        
        [x,y] = ginput;
        
        for k = 1:numel(x)
            %create 15x15 box around that point
            crop = img(round(y(k))-7:round(y(k))+7, round(x(k))-7:round(x(k))+7);
            %To find background, take pixels 4 deep around edge
            %Take avg of bottom 25 percent
            seg1 = crop(1:4,:);
            seg2 = crop(12:15,:);
            seg3 = crop(5:11,1:4);
            seg4 = crop(5:11,12:15);
            s = sort([seg1(:) ; seg2(:) ; seg3(:) ; seg4(:)]);
            bkgrd = mean(s(1:44))
            max(crop(:))
            intensity = max(crop(:)) - bkgrd;
            intData(m,:) = [i j k intensity];
            m = m+1;
        end
        clf; pause(.1);
    end
end

dlmwrite('intron_intensity.txt',intData,'\t');