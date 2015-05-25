load('/Volumes/Data2_OPM/OrganizedData/A549/130611_TBCB/data017.mat' ;

obj = objects(1);
channel = 'nir';

coords = obj.channels.(channel).spotCoordinates;
x = coords(:,2);
y = coords(:,1);
z = coords(:,3);

dapiMask = obj.channels.dapi.processor.mask;

tmp1 = max(dapiMask,[],1); % same dim as x
tmp2 = max(dapiMask,[],2); % same dim as y

xmin = min(find(tmp1>0));
xmax = max(find(tmp1>0));
xmean = (xmin+xmax)/2;

ymin = min(find(tmp2>0));
ymax = max(find(tmp2>0));
ymean = (ymin+ymax)/2;

window = 10;

xsub = coords((y>(ymin-window) & y<(ymin+window)),:);
xsubx = xsub(:,2);
xsubz = xsub(:,3);
ysub = coords((x>(xmin-window) & x<(xmin+window)),:);
ysuby = ysub(:,1);
ysubz = ysub(:,3);

%cellOutline_140115;

outline_y = celltop(round(ymean),:);
outline_x = celltop(:,round(xmean));

%outline_ymax = max(celltop,[],1);
%outline_xmax = max(celltop,[],2);

figure; scatter(xsubx,xsubz); hold on; line([xmin xmin],[min(xsubz) max(xsubz)]); line([xmax xmax],[min(xsubz) max(xsubz)]); plot(1:numel(outline_y),outline_y,'r');
figure; scatter(ysuby,ysubz); hold on; line([ymin ymin],[min(ysubz) max(ysubz)]); line([ymax ymax],[min(ysubz) max(ysubz)]); plot(1:numel(outline_x),outline_x,'r');

%figure; scatter(x,z); hold on; line([xmin xmin],[min(z) max(z)]); line([xmax xmax],[min(z) max(z)]);
%figure; scatter(y,z); hold on; line([ymin ymin],[min(z) max(z)]); line([ymax ymax],[min(z) max(z)]);

