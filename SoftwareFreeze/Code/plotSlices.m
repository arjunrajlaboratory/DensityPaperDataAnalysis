%% 

figure; imshow(objmask); hold on;
line([1 size(objmask,2)],[500 500],'LineWidth',2,'Color','r','LineSmoothing','on');
%line([1 size(objmask,2)],[200 200],'LineWidth',2,'Color','r','LineSmoothing','on');
%line([1 size(objmask,2)],[300 300],'LineWidth',2,'Color','r','LineSmoothing','on');
%line([1 size(objmask,2)],[400 400],'LineWidth',2,'Color','r','LineSmoothing','on');
%line([1 size(objmask,2)],[500 500],'LineWidth',2,'Color','r','LineSmoothing','on');

%figure; imshow(objmask); hold on;
%line([40 40],[1 size(objmask,1)],'LineWidth',2,'Color','r','LineSmoothing','on');
%line([80 80],[1 size(objmask,1)],'LineWidth',2,'Color','r','LineSmoothing','on');
%line([120 120],[1 size(objmask,1)],'LineWidth',2,'Color','r','LineSmoothing','on');
%line([160 160],[1 size(objmask,1)],'LineWidth',2,'Color','r','LineSmoothing','on');
%line([200 200],[1 size(objmask,1)],'LineWidth',2,'Color','r','LineSmoothing','on');

%% This takes slices moving down through y

idx = 200;
num = 1;

topline = celltop(idx,1:size(celltop,2));
sidxtop = min(find(topline~=0));
eidxtop = max(find(topline~=0));
topline(find(topline==0)) = [];
figure; 
line(sidxtop:eidxtop, topline,'LineWidth',2,'Color','k','LineSmoothing','on'); hold on;
%line(1:size(celltop,2),topline); hold on;

bottomline = cellbottom(idx,1:size(cellbottom,2));
sidxbottom = min(find(bottomline~=0));
eidxbottom = max(find(bottomline~=0));
bottomline(find(bottomline==0)) = [];
line(sidxbottom:eidxbottom, bottomline,'LineWidth',2,'Color','k','LineSmoothing','on');

line(1:size(objmask,2),objmask(idx,1:size(objmask,2))+9,'LineWidth',1,'LineSmoothing','on');

%Draws line denoting area blocked by mask
% msk = objmask(idx,1:size(objmask,2));
% d = diff(msk);
% didx = find(d~=0);
% 
% didx = didx(didx > sidxtop + 1);
% didx = didx(didx < eidxtop - 1);
% 
% if(didx)
%     for i = 1:numel(didx)
%         line([didx(i) didx(i)],[topline(didx(i)-didx(1)+1) bottomline(didx(i)-didx(1)+1)],'LineWidth',2,'Color','k');
%     end
% end

%to connect
line([sidxtop sidxbottom],[topline(1) bottomline(1)],'LineWidth',2,'Color','k','LineSmoothing','on');
line([eidxtop eidxbottom],[topline(end) bottomline(end)],'LineWidth',2,'Color','k','LineSmoothing','on');

% idxlist = horzcat(find(yp>idx+num), find(yp<idx-num));
% xlist = xp;
% xlist(idxlist) = [];
% zlist = zp;
% zlist(idxlist) = [];

% scatter(xlist,zlist,'filled','m');

idxlist = horzcat(find(xp>idx+num), find(xp<idx-num));
ylist = yp;
ylist(idxlist) = [];
zlist = zp;
zlist(idxlist) = [];

scatter(ylist,zlist,'filled','m');

%axis([0 size(objmask,2) 8.5 max(zp)+.5]);
%axis off;

%% This takes slices moving left to right through x

idx = 400;
num = 5;

topline = celltop(1:size(celltop,1),idx);
sidxtop = min(find(topline~=0));
eidxtop = max(find(topline~=0));
topline(find(topline==0)) = [];
figure; 
line(sidxtop:eidxtop, topline,'LineWidth',2,'Color','k','LineSmoothing','on'); hold on;

bottomline = cellbottom(1:size(cellbottom,1),idx);
sidxbottom = min(find(bottomline~=0));
eidxbottom = max(find(bottomline~=0));
bottomline(find(bottomline==0)) = [];
line(sidxbottom:eidxbottom, bottomline,'LineWidth',2,'Color','k','LineSmoothing','on');

%to connect
line([sidxtop sidxbottom],[topline(1) bottomline(1)],'LineWidth',2,'Color','k','LineSmoothing','on');
line([eidxtop eidxbottom],[topline(end) bottomline(end)],'LineWidth',2,'Color','k','LineSmoothing','on');

line(1:size(objmask,1),objmask(1:size(objmask,1),idx)+8,'LineWidth',1,'LineSmoothing','on');

idxlist = horzcat(find(xp>idx+num), find(xp<idx-num));
ylist = yp;
ylist(idxlist) = [];
zlist = zp;
zlist(idxlist) = [];

scatter(ylist,zlist,'filled','m');

axis([0 size(objmask,1) 7.5 max(zp)+.5]);
axis off;