contents = dir('data*');

x = load('volume_expand_battleship.txt');

for i = 1:size(x,1)
    load(contents(x(i,1)).name);
    objects(x(i,2)).metadata.volume = x(i,3);
    save(contents(x(i,1)).name,'objects');
end

%%

% contents = dir('data*');
% color = 'alexa';
% nSpots = 0;
% m = 0;
% 
% for i=1:numel(contents)
%     load(contents(i).name);
%     for j=1:numel(objects)
%         nSpots = nSpots + numel(objects(j).channels.(color).spotCoordinates(:,1));
%         m = m+1;
%     end
% end
% 
% avg = nSpots/m