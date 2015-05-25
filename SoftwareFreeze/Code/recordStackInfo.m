function [] = recordStackInfo(varargin)

%Records the plane spacing within stacks to objects().metadata.planeSpacing
%recordStackInfo(spacing) applies the same spacing to all objects
%recordStackInfo(spacing,objs) applies spacing to specified objects

contents = dir('data*');

spacing = varargin{1};

if nargin > 1
    objs = varargin{2};
else
    objs = 1:numel(contents);
end

for i = objs
    disp(['loading' contents(i).name]);
    load(contents(i).name);
    for j = 1:numel(objects)
        objects(j).metadata.planeSpacing = spacing;
    end
    disp(['saving' contents(i).name]);
    save(contents(i).name,'objects')
end