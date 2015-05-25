dayContents = dir('13*beads');
mainDir = pwd;
areaTable = 0;
heightTable = 0;
m = 0;

for i = 1:numel(dayContents) % load data taken on different days
    cd([mainDir '/' dayContents(i).name]);
    wellContents = dir('well*');
    if numel(wellContents) > 3
        range = [6 5 4 3 2 1];
    else
        range = [3 2 1];
    end
    wellDir = pwd;
    
    m = 0;
    
    for j = range % load conditions: live, fixed, ethanol
        
        if mod(j,3) == 0
            if areaTable == 0
                startval = 0;
            else
                startval = size(areaTable,2);
            end
        end
        
        m = m+1;
        cd([wellDir '/' wellContents(j).name]);
        %disp([wellDir '/' wellContents(j).name]);
        dataDir = dir('data*');
        
        for k = 1:numel(dataDir)
            load(dataDir(k).name);
            areaTable(mod(m-1,3)+1,k+startval) = numel(find(objects(1).object_mask.mask));
            heightTable(mod(m-1,3)+1,k+startval) = objects(1).metadata.cellTop-objects(1).metadata.cellBottom;
            %disp(k)
        end
    end
end

figure; boxplot(heightTable','labels',{'live', 'fixed', 'ethanol'}); title('heights');
%figure; boxplot(areaTable'); title('areas');
figure; scatter(areaTable(1,:),areaTable(2,:)); title('areas'); xlabel('live'); ylabel('fixed');
figure; scatter(areaTable(1,:),areaTable(3,:)); title('areas'); xlabel('live'); ylabel('ethanol');
figure; scatter(areaTable(2,:),areaTable(3,:)); title('areas'); xlabel('fixed'); ylabel('ethanol');
figure; bar(heightTable'); title('heights');

%%%%%%

% c = dir('well2*');
% wd = pwd;
% k = 0;
% 
% areaTable = 0;%zeros(numel(cs),3);
% heightTable = 0;%zeros(numel(cs),3);
% 
% %       live fixed ethanol
% % cell1
% % cell2
% %.
% %.
% %.
% 
% for i = [3 2 1]
%     k = k+1;
%     cd([wd '/' c(i).name]);
%     cs = dir('data*');
% 
%     for j = 1:numel(cs)
%         load(cs(j).name);
%         areaTable(k,j) = numel(find(objects(1).object_mask.mask));
%         heightTable(k,j) = objects(1).metadata.cellTop-objects(1).metadata.cellBottom;
%     end
% end
