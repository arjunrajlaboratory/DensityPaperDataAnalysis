plotrows = 2;
plotcols = 3;


figure(3);clf;
subplot(plotrows,plotcols,1)
%imshow(obj.channels.(rna_color).maxlaplaceimage)
%subplot(1,1,1)
imshow(obj.channels.alexa.maxlaplaceimage)
hold on

heightmap=jet(max(zp)-min(zp)+1);

for i=1:length(zp)

line(xp_adjusted(i),yp_adjusted(i),'Marker','.','MarkerEdgeColor',...
    heightmap(zp(i)-min(zp)+1,:))
end

hold off

figure(3);subplot(plotrows,plotcols,2)
imshow(obj.channels.alexa.maxlaplaceimage)
hold on
dt=DelaunayTri(xp_tent',yp_tent');
hold on
triplot(dt)    

for i=1:length(zp_tent)
line(xp_tent(i),yp_tent(i),'Marker','.','MarkerEdgeColor',...
    heightmap(zp_tent(i)-min(zp)+1,:))
end
for i=1:length(zp_out)
line(xp_out(i),yp_out(i),'Marker','.','MarkerEdgeColor',...
    heightmap(zp_out(i)-min(zp)+1,:))
end
%plot(xp_tent(hullindices),yp_tent(hullindices),'r-')
hold off

figure(3);subplot(plotrows,plotcols,3)
imshow(obj.channels.alexa.maxlaplaceimage)
hold on
dt=DelaunayTri(xp_tent',yp_tent');
hold on
triplot(dt)

for i=1:length(zp_tent)
line(xp_tent(i),yp_tent(i),'Marker','.','MarkerEdgeColor',...
    heightmap(zp_tent(i)-min(zp)+1,:))
end
for i=1:length(zp_out)
line(xp_out(i),yp_out(i),'Marker','.','MarkerEdgeColor',...
    heightmap(zp_out(i)-min(zp)+1,:))
end


% Need to now create topography maps. We can start with plotting the mask

figure(3);subplot(plotrows,plotcols,4)

imshow(obj.object_mask.mask)

% ideally for every pixel I would be able to calculate a height.

% It looks like TriScatteredInterp is just the tool for the job.



%% TOpography plot


figure(4);clf;subplot(1,plotcols,1)
imshow((topZs-min(zp))/max(zp));colormap jet
subplot(1,plotcols,2)
imshow((bottomZs-min(zp))/max(zp));colormap jet

subplot(1,plotcols,3)
heightZs=topZs-bottomZs;
imshow((heightZs)/(max(zp)-min(zp)));colormap jet


heightZs(~(obj.object_mask.mask))=0;
heightZs(isnan(heightZs))=0;


%% Topography plot again


%figure(4);clf;subplot(1,2,1)
%imshow((topZs-min(zp))/max(zp));colormap jet;

%Ftopalt=TriScatteredInterp(toptent(:,1),toptent(:,2),toptent(:,3),'natural');

%subplot(1,2,2)
%topZsalt=Ftopalt(Xs,Ys);
%imshow((topZsalt-min(zp))/max(zp));colormap jet;
