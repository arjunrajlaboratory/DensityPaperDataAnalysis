intColor = 'tmr';
rnaColor = 'alexa';
contents = dir('data*');
m = 1;
n = 1;
clear intData;
clear avgData;

for i = 1:numel(contents)
    load(contents(i).name);
    fprintf('loading %s\n',contents(i).name);
    for j = 1:numel(objects)
        
        obj = objects(j);
        avgInt = 0;
        
        %Returns indices of tmr spots in nucleus
        dapiMask = obj.channels.dapi.processor.mask;
        a1 = obj.channels.(intColor).inNucleus(dapiMask);
        
        %Returns indices of tmr spots colocalized with alexa spots
        a2 = obj.channels.(intColor).getColocalized(obj.channels.(rnaColor),3,.35/.125);
        
        %Intron spot coordinates
        coords = obj.channels.(intColor).spotCoordinates(a1&a2,:);
        rc = round(coords);
        
        %coords = obj.channels.(intColor).spotCoordinates;
        bbox = obj.object_mask.boundingbox; %upperleft_x,upperleft_y, x_width, y_width
        
        %Keep only spots that are in the nucleus
        %rc = round(coords);
        %x = rc(:,2);
        %y = rc(:,1);
        %ind = sub2ind(size(dapiMask),y,x);
        
        %rc(find(dapiMask(ind)==0),:) = [];
        
        for k=1:size(rc,1)
            %Locate the correct plane & crop to get same region as mask
            stk = obj.channelStk(intColor);
            img = stk(:,:,rc(k,3));
            
            %Crop the image to only a 15x15 box around spot
            crop = img(rc(k,1)-7:rc(k,1)+7, rc(k,2)-7:rc(k,2)+7);
            
            %To find background, take pixels 4 deep around edge
            %Take avg of bottom 25 percent
            seg1 = crop(1:4,:);
            seg2 = crop(12:15,:);
            seg3 = crop(5:11,1:4);
            seg4 = crop(5:11,12:15);
            s = sort([seg1(:) ; seg2(:) ; seg3(:) ; seg4(:)]);
            bkgrd = mean(s(1:44));
            intensity = max(crop(:)) - bkgrd;
            intData(m,:) = [i j k intensity];
            avgInt = avgInt + intensity;
            m = m+1;
        end
        %avgData(n,:) = [i j avgInt/size(rc,1)];
        avgData(n,:) = [i j size(rc,1) avgInt];
        n = n+1;
    end
end

dlmwrite('intron_intensity.txt',intData,'\t');
dlmwrite('intron_intensity_avg.txt',avgData,'\t');