% Locate z position for each transcription site, using intron channel

intronChannel = 'tmr';

contents = dir('data*');

m=1;

for i = 1:numel(contents) % each data file
    load(contents(i).name);
    for j = 1:numel(objects) % each image object
        obj = objects(j);
        if ~isfield(obj.metadata,'transcriptionSitesManual')
            continue;
        end
        im = obj.channelStk(intronChannel);
        Xs = obj.metadata.transcriptionSitesManual.Xs;
        Ys = obj.metadata.transcriptionSitesManual.Ys;
        xCoords = round(Xs);
        yCoords = round(Ys);
        for k = 1:numel(obj.metadata.transcriptionSitesManual.Xs) % each transcription site
            x = xCoords(k);
            y = yCoords(k);
            imCropTiny = im(y,x,:);
            [val plane] = max(imCropTiny);
            imCropLarge = im((y-10):(y+10),(x-10):(x+10),plane);
            
            % new
            imLarge = im(:,:,plane);
            imLargeFilt = medfilt2(imLarge);
            regMax = imregionalmax(imLargeFilt);
            indices = find(regMax);
            [coordsX coordsY] = ind2sub(size(regMax),indices);
            idx = knnsearch([coordsX coordsY],[y x]);
            maxVal = imLarge(coordsX(idx), coordsY(idx));
            intensity = maxVal - min(imCropLarge(:));
            %intensity = max(imCropLarge(:)) - min(imCropLarge(:));
            % end new
            
            data(m,:) = [i j intensity];
            m = m + 1;
        end
    end
end


load('data010.mat')
obj = objects(1);

Xs = obj.metadata.transcriptionSitesManual.Xs;
Ys = obj.metadata.transcriptionSitesManual.Ys;

xCoords = round(Xs);
yCoords = round(Ys);

im = obj.channelStk('tmr');
imCropTiny = im(yCoords(1),xCoords(1),:);
[val plane] = max(imCropTiny);

imLarge = im(:,:,plane);
imLargeFilt = medfilt2(imLarge);
regMax = imregionalmax(imLargeFilt);

indices = find(regMax);
[coordsX coordsY] = ind2sub(size(regMax),indices);

idx = knnsearch([coordsX coordsY],[y x]);

imCropLarge = im((yCoords(1)-10):(yCoords(1)+10),(xCoords(1)-10):(xCoords(1)+10),plane);

maxVal = imLarge(coordsX(idx), coordsY(idx));
%maxVal = max(imCropLarge(:));
minVal = min(imCropLarge(:));

intensity = maxVal - minVal;