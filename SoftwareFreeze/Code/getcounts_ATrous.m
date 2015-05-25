%%
contents = dir('data*');

fillChannel = 'nir';
countChannel = 'alexa';
cyclinChannel = 'cy';

for i = 1:numel(contents)
    fprintf('Loading %s\n',contents(i).name);
    load(contents(i).name);
    for j = 1:numel(objects)
        
        %for ACTN4
        %if i==4 & j==1
        %    continue;
        %elseif i==6 & j==2
        %    continue;
        %end
        
        %if i==10 & j==1
        %    continue;
        %end

        
        %Find spot locations of 'fill' channel
        [x y z] = findSpotsATrousCC(objects(j),fillChannel,'lenient');
        objects(j).metadata.(fillChannel).x = x;
        objects(j).metadata.(fillChannel).y = y;
        objects(j).metadata.(fillChannel).z = z;
        %n1 = numel(x);
        
        %[x y z] = findSpotsATrousCC(objects(j),cyclinChannel,'normal');
        %objects(j).metadata.(cyclinChannel).x = x;
        %objects(j).metadata.(cyclinChannel).y = y;
        %objects(j).metadata.(cyclinChannel).z = z;

        %[x y z] = findSpotsATrousCC(objects(j),countChannel,'lenient');
        %objects(j).metadata.(countChannel).x = x;
        %objects(j).metadata.(countChannel).y = y;
        %objects(j).metadata.(countChannel).z = z;
        %n2 = numel(x);
        
    end
    fprintf('Saving %s\n',contents(i).name);
    save(contents(i).name,'objects');
end

