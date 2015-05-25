%%
contents = dir('data*');
m = 1;
numRNA = 0;
color = 'alexa';

for i = 1:numel(contents)
    load(contents(i).name);
    for j = 1:numel(objects)
        numRNA = numRNA + objects(j).channels.(color).numSpots;
        m = m+1;
    end
end

avgRNA = numRNA/m