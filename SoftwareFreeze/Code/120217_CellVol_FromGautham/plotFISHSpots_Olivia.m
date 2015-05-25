%% plotFISHSpots  
% plots the 3D coordinates of FISH spots stored in |image_object|.
% This function can plot points for both a single |image_object| in it's own
% BoundingBox, or an array of |image_object| all within a microscope 
% field-of-view (FOV) defined as a 1024 x 1024 x zMax 3D grid. 
%

%% Input
% *|image_object|* - the @image_object instance
%
% *|zMax|* - this is the height of the stack (number of planes) used 
% to draw plot (double)
%
% _Optional_
%
% *|excludeFields|* - cell array of strings listing the channel names to be
% exlcuded from being plotted. (eg. {'tmr','alexa'})
%
% *|chrShiftLib|* - structure with fields containing interpolant functions that 
% map XYZ coordinates to correct for chromatic shift.
%
%  (eg: chrShiftLib -> alexa (to nir) -> dX() & dY() & dZ()
%                      gfp   (to nir) -> dX() & dY() & dZ()                
%                      tmr   (to nir) -> dX() & dY() & dZ()
%
% *|DAPIFlag|* - when using a single |image_object|, DAPI outline is found 
% and plotted in the XY plane for 2D viewing (BOOLEAN)
%
% *|scatterAxes|* - plot in a provided axes handle instead of in a new figure 
%

%% Example Usage:
%  >> load data001; % load the image_objects stored in current working directory
%  >> plotFISHSpots(objects,zMax,'excludeFields',{'alexa'});

%% Notes
% There is a *SPECIAL FEATURE* in this function: If the excludeFields
% specifies channel names to be excluded and it results in only 2 remaining
% fields, the CHROMATIC SHIFT between multi-colored points is calculated 
% and visualized in the plot.
% 
% This function uses |scatter3| since |plot3| cannot do points in more than one
% color successfully (try it yourself if you don't believe me).
% A new figure window is open to draw the 3D scatter plot.


function [varargout] = plotFISHSpots(imgObjs,zMax,varargin)

 
%----------------------------------------------
% Parse input parameters and load the STK files
%----------------------------------------------

p = inputParser;
p.addOptional('spotSource','fitdataRNAonly',...
                    @(x) any(strcmpi(x,{'fitdataRNAonly','fitdata','pseudoColorSpot'})));
p.addOptional('excludeFields',{},@iscell);  % cell array of strings
p.addOptional('chrShiftLib',struct(),@isstruct); % struct w/ chromatic shift maps
p.addOptional('DAPIFlag',false,@islogical); % plot the outline of DAPI on the Z=0 plane
p.addOptional('scatterAxes',[],@ishandle);
p.addOptional('markerStyle','skittles',@(x) any(strcmpi(x,{'skittles','fruitloops',...
                                                           'pluses','triangles'})));
p.addOptional('viewingAngle','2D',@(x) any(strcmpi(x,{'2D','3D'})));
p.parse(varargin{:});

spotSource    = p.Results.spotSource;
excludeFields = p.Results.excludeFields;
chrShiftLib   = p.Results.chrShiftLib;
DAPIFlag      = p.Results.DAPIFlag;
scatterAxes   = p.Results.scatterAxes;
markerStyle   = p.Results.markerStyle;
viewingAngle  = p.Results.viewingAngle;

% If we are working with an entire field of imgObjs (aka: an array of 
% imgObjs loaded from data00*.mat) Use the first one as a reference to
% get channels used for introns.
xyLimits = [0 1024 0 1024];
if length(imgObjs) == 1   
    imgObj = imgObjs;
    % With only one image_object, plot w/in BoundingBox. Axis coordinates is IJ so origin is top left
    xyLimits = imgObj.object_mask.boundingbox;  % [ upper-left corner X,Y , x_width, y_width]
    % *NOTE* : x-y are swapped in segmenting
    xyLimits = [xyLimits(1) xyLimits(1)+xyLimits(3) ...
                xyLimits(2) xyLimits(2)+xyLimits(4)]; % [xmin xmax ymin ymax]
else
    imgObj = imgObjs(1);
end

plotChannels = imgObj.RNAFields;  % cell array of strings
% Remove the fields that do not contain spot coordinate data. Also remove fields
% specified in 'excludeFields' input variable
for chan = imgObj.RNAFields  
    chan = char(chan);  % convert cell -> str
    if ~isempty(strfind(spotSource,'fitdata')) ...
             && ~isfield(imgObj.channels.(chan),spotSource) % no fitdata 
        plotChannels = plotChannels(~strcmp(chan,plotChannels));  % remove channel
    end

    if any(strcmp(chan,excludeFields))  % remove explicitly stated channels
        plotChannels = plotChannels(~strcmp(chan,plotChannels));  % remove chan
    end
end

if strcmp(spotSource,'pseudoColorSpot')
    plotChannels = [plotChannels {'ref'}];
end
%^^^^^^^^^^^^^^^^^ Done with initial loading of data ^^^^^^^^^^^^^^^^^^^^^


%-------------------------------------------------------------------------
% Translate FISH spots from their image_object coordinates into the full
% microscope field coordinates. Also handle adjustment of chromatic shift
% if a reference library was defined for each channel.
%-------------------------------------------------------------------------

% We will construct [m x 3] XYZ 3D coordinate data for each intron channels 
plotPoints = struct(); 
for chan = plotChannels     % for each intron channel
    chan = char(chan);        % convert cell -> str
    plotPoints.(chan) = []; % instantiate field with empty array
end

for imgObj = imgObjs  % for each image_object in the input array imgObjs

    if isempty(scatterAxes)
        % spot coordinates are stored in reference to the bounding box of the
        % image_object mask defined during segmentation. We need to convert to 
        % microscope field coordinates. * MAKE SURE TO WATCH THE X-Y SWAP! *  Axis IJ
        box = imgObj.object_mask.boundingbox;  % [ upper-left corner X,Y , x_width, y_width] 
        plusX = box(2); plusY = box(1);   % XY swap is bounding box IJ = XY vs img_proc IJ = YX
        %plusX = box(2); plusY = box(1);  % *NOTE* : x-y are swapped in segmenting
    else
        % When sending plot to a provided |axes| handle, keep the spots in their local 
        % |image_object| coordinates
        plusX = 0; plusY = 0; 
    end
    
    if strfind(spotSource,'fitdata')
        for chan = plotChannels % for each channel in image_object containing spot data
            chan = char(chan);  % convert cell -> str
            
            chrShift = []; % Check if we have a shift library for this img channel
            if isfield(chrShiftLib,chan); chrShift = chrShiftLib.(chan); end;

            for i = 1:length(imgObj.channels.(chan).(spotSource).xp_fit) % for each spot
                xyz = [imgObj.channels.(chan).(spotSource).xp_fit(i)...  % imgObj coordinates
                       imgObj.channels.(chan).(spotSource).yp_fit(i)...  
                       imgObj.channels.(chan).(spotSource).rawzp(i) ]; % return array [x y z]
                
                if ~isempty(chrShift)  % If we have a mapping fxn
                    xyz = [xyz(1)+chrShift.dX(double(xyz)) ...  % correct the
                           xyz(2)+chrShift.dY(double(xyz)) ...  % chromatic
                           xyz(3)+chrShift.dZ(double(xyz))];    % shift !
                end

                xyz = xyz + [plusX plusY 0];  % convert to FOV coordinates

                plotPoints.(chan) = [plotPoints.(chan); xyz];
            end
        end

    elseif strcmp(spotSource,'pseudoColorSpot')
        for pSpot = imgObj.pseudoColorSpots
            if pSpot.isShiftReference
                chan = 'ref';
            else
                randomColorInd = ceil(rand*(length(plotChannels)-1)); % exclude the reference index
                chan = char(plotChannels(randomColorInd));
            end
            xyz = pSpot.position + [plusX plusY 0];
            plotPoints.(chan) = [plotPoints.(chan); xyz]; 
        end

    end
end

%----------------
% PLOT the points
%----------------
shift = [];
if isempty(fieldnames(plotPoints))
    fprintf(1,'No channels had FISH spot data, no points to plot\n');
    return
end

% When a scatterAxes is provided, we do not draw to a new figure
if isempty(scatterAxes)
    figure;
    scatterAxes = axes;
end
scatterH = []; % Hold all the ScatterGroup object handles
hold(scatterAxes,'on')  % We iteratively add the ScatterGroup objects for each channel to the figure
colors.alexa = [1 0 1]; % magenta
colors.tmr   = [0 1 0]; % Green
colors.gfp   = [0 1 1]; % Cyan
colors.nir   = [1 0 0]; % Red
colors.cy    = [0 0 1]; % Blue
colors.ref   = [0 0 0]; % Black for spots used in chromatic shift correction
legendTitles = {};
for chan = plotChannels % for each channel, plot the spots using scatter3()
    chan = char(chan);  % convert cell -> str
    XYZ = plotPoints.(chan);
    if isempty(XYZ); continue; end;
    C = [];for i = 1:size(XYZ,1); C = [C;colors.(chan)]; end; % construct color array
    % Swap XY so that our spot coordinates IJ = YX is converted to IJ = XY
    H = scatter3(XYZ(:,2),XYZ(:,1),XYZ(:,3),50,C,'Parent',scatterAxes); 
    if strcmpi(markerStyle,'skittles')
        set(H,'MarkerFaceColor',colors.(chan),'MarkerEdgeColor','black');  % fill the circles
    elseif strcmpi(markerStyle,'fruitloops')
        set(H,'LineWidth',2,'SizeData',100);  % make big colored rings
    elseif strcmpi(markerStyle,'pluses')
        set(H,'Marker','+','LineWidth',2,'SizeData',100);
    elseif strcmpi(markerStyle,'pyramid')
        set(H,'Marker','^','LineWidth',2,'SizeData',100);
    end
    scatterH     = [scatterH; H];  % store the handles of the ScatterGroups
    legendTitles = [legendTitles; {chan}];  % store the channel names for the legend
end
% Formatting and legend labels
grid(scatterAxes,'on');

if DAPIFlag
    if ~isfield(imgObj.metadata,'dapiBorder')
        dapiStk = imgObj.channelStk('dapi');
        dapiStk = max(dapiStk,[],3);
       % dapiStk(~imgObj.object_mask.mask) = min(dapiStk(:));
        dapiMask = maskWithDapi(dapiStk);
        B = bwboundaries(dapiMask);
    else
        B = imgObj.metadata.dapiBorder;
    end
    plusXB = repmat(plusX,length(B{1}),1);
    plusYB = repmat(plusY,length(B{1}),1);
    DAPIlineH = plot(B{1}(:,2)+plusYB,B{1}(:,1)+plusXB,'r','Parent',scatterAxes);  
end
xlim([xyLimits(1) xyLimits(2)]); xlabel('x','FontSize',16);
ylim([xyLimits(3) xyLimits(4)]); ylabel('y','FontSize',16);
zlim([0 zMax]); zlabel('z','FontSize',16);
set(scatterAxes,'YDir','reverse');
legend(scatterH,legendTitles{:});

if strcmpi(viewingAngle,'3D')
    view(scatterAxes,45,20);  % Tilt the view so it looks more "Three-Dee"
    set(scatterAxes,'Projection','perspective');
elseif strcmpi(viewingAngle,'2D')
    view(scatterAxes,2);  % top-down 2D view
end

    

% *SPECIAL FEATURE* - see heading for info
if length(plotChannels) == 2 && ~isempty(excludeFields) 
% Find the vectors representing the shift in X-Y-Z between point sets A & B
    fromPts = double(plotPoints.(char(plotChannels(2)))); 
    toPts   = double(plotPoints.(char(plotChannels(1))));
    shift = toPts - fromPts;
    % draw velocity field arrows pointing from points A to B, no scaling
    quiver3(fromPts(:,2),fromPts(:,1),fromPts(:,3),...  % XYZ positions
                  shift(:,2),shift(:,1),shift(:,3),...  % UVW vectors
                  'Parent','scatterAxes',0,'LineWidth',2);  % 1X scaling of magnitude 
end
hold(scatterAxes,'off');

if nargout == 1
    varargout = {scatterH};
elseif nargout == 2
    if ~DAPIFlag
        fprintf(1,'DAPIFlag was not set so no DAPI plot was returned');    
        DAPILineH = [];
    end
    varargout = [{scatterH},{DAPIlineH}];
end
% ALL DONE
