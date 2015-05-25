wd = pwd;

contents = dir('well*');

for i = 1:numel(contents)
    cd([wd '/' contents(i).name]);
    datafolders = dir('data*');
    if i == 1
        areas = zeros([numel(datafolders) 3]);
    end
    for j = 1:numel(datafolders)
        load(datafolders(j).name);
        obj = objects(1);
        areas(j,i) = numel(find(obj.object_mask.mask));
    end
end

dlmwrite('~/Dropbox/densitypaper/SupplementaryFigures/FixationAndVolume/Data_130115_2.txt',areas,'delimiter','\t','precision',8);