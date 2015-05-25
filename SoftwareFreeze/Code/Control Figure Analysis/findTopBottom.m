c = dir('well*');
wd = pwd;

for i = 1:numel(c)
    %text = sprintf('Changing directory to %c\n',c(i).name);
    disp(['Changing directory to ' c(i).name]);
    cd([wd '/' c(i).name])
    cs = dir('data*');
    for j = 1:numel(cs)
        load(cs(j).name)
        for k = 1:numel(objects)
            if isfield(objects(k).metadata,'cellBottom')
                continue;
            end
            StackViewer(objects(k).channelStk('cy'))
            bottom = input('plane number of bottom of cell: ');
            top = input('plane number of top of cell: ');
            close StackViewer;
            objects(k).metadata.cellBottom = bottom;
            objects(k).metadata.cellTop = top;
        end
        save(cs(j).name,'objects');
    end
end