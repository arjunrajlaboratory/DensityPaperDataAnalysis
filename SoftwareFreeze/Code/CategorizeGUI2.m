%% SegmentGUI
% A graphical user interface to segment image objects (cells, embryos, nuclei, etc)

%% Description
% This program lets you select the cells or embryos for which you would like
% to count RNA spots.
%
% Press the "Segment" button or hit the 'option' key to use the freehand
% segmenting tool.
%
% Press the "Undo Segment" button or hit the 'delete' key to undo the last
% segmented object.
%
% Press the "Next File" button or hit the 'rightarrow' key to save and
% move on to the next file. "Prev File" ('left arrow') to save and go
% back in the file sequence.
%

%% MATLAB Figure file for GUI development
% There is a MATLAB figure file (SegmentGUILayout.fig) that is only used as
% a way to place |uicontrol()| elements and use position properties field
% values to layout the GUI programmatically in this m-file.

%% Filename Enforcement
% Uniform file name schemes are used in the Raj lab for image stack files.
% Fluorescent channel images use short names followed by 3-digit number.
% The used names are {'tmr','alexa','cy','gfp','nir'}
% DAPI (nuclear stain) and transmission images should follow the same
% numbering scheme.

%% Credits
% * Originally conceived by Arjun Raj
% * GUI-fied by Yaanik Desai
% * Rewritten by Marshall J. Levesque 2012

function varargout = CategorizeGUI2(varargin)

%======================================================
% Global variables are all stored in "Hs" structure,
% aka "handles", the common place to store app data
% in a MATLAB GUI program. "Hs" makes for cleaner code
% and "handles" doesn't apply to all the data types we
% store in it.
%======================================================

figH = figure('Position',[200 100 806 567],...
    'NumberTitle','off',...
    'Name','CategorizeGUI',...
    'Resize','on',...
    'Toolbar','none',...
    'MenuBar','none',...
    'Color',[0.247 0.247 0.247],...
    'KeyPressFcn',@KeyPressCallback,...
    'Visible','off');

Hs = guihandles(figH);
Hs.figH = figH;

% Initialize some bookkeeping variables
Hs.progressCount = 0;

% Start by checking if we already have 'data***.mat' files
% that store |image_object|. If so, we load the data files
% to get images and the segmentation ROIs.
%---------------------------------------------------------
Hs.dirPath = pwd;
[Hs.dataFiles,Hs.dataNums] = getDataFiles(Hs.dirPath);

% Look for image files in the current working directory. If we
% don't find any, ask the user to navigate to the image files
% using the GUI file browser. We check again for 'data***.mat'
% and then the image files
%--------------------------------------------------------------
if isempty(Hs.dataFiles)  % no image files here
    fprintf(1,'Could not find any image files!\n');
    yn = input('Navigate to the directory with your files? y/n [y]','s');
    if isempty(yn); yn = 'y'; end;  % default answer when press return only
    if any(strcmp(yn,{'y','Y','yes','Yes','YES','1'}))
        Hs.dirPath = uigetdir(pwd,'Navigate to image files');
        if Hs.dirPath == 0;  % User pressed cancel
            return;  % quit the GUI
        end
        [Hs.dataFiles,Hs.dataNums] = getDataFiles(Hs.dirPath);
        [Hs.foundChannels,Hs.fileNums,Hs.imgExts] = getImageFiles(Hs.dirPath);
        if isempty(Hs.fileNums)
            error('Could not find image files (eg ''tmr001.tif'' etc) to segment');
        end
    else
        return;  % user did not want to navigate to image files, quit GUI
    end
end

Hs.currentFile = 1;
Hs.currentObject = 1;

set(figH,'Visible','on');
Hs = buildGUI(Hs,figH);
Hs.imgH = [];
Hs = makeAndShowOverlay(Hs);

% Update handles structure
%--------------------------
guidata(figH, Hs);

return;
%============================================================
% End of the main SegmentGUI function. Everything else happens
% passing around the |Hs| structure between functions
%============================================================


function CloseRequest(hObject,eventdata)
Hs = guidata(gcbo);
setappdata(0,'Success',true);
guidata(gcbo,Hs);
uiresume(hObject);
delete(hObject);




function Hs = makeAndShowOverlay(Hs)

load(Hs.dataFiles(Hs.currentFile).name);
alexaImage = objects(Hs.currentObject).channels.alexa.processor.getImage;
alexaRGB(:,:,1) = alexaImage;
alexaRGB(:,:,2) = alexaImage;
alexaRGB(:,:,3) = alexaImage;

dapiImage = objects(Hs.currentObject).channels.dapi.processor.getImage;
dapiRGB(:,:,1) = .556*dapiImage;
dapiRGB(:,:,2) = .2078*dapiImage;
dapiRGB(:,:,3) = .93725*dapiImage;


object_mask = objects(Hs.currentObject).object_mask.mask;
object_mask = bwperim(object_mask);
object_mask = imdilate(object_mask, strel('line', 1, 0));
object_mask = imdilate(object_mask, strel('line', 1, 90));
object_mask = single(object_mask);
color_mask(:,:,1) = object_mask;
color_mask(:,:,2) = object_mask;
color_mask(:,:,3) = object_mask;

Hs.RGB = scale(alexaRGB + dapiRGB + color_mask);

Hs.imgH = imshow(Hs.RGB,'Parent',Hs.imgAx);
set(Hs.imgAx,'XLim',[1 size(Hs.RGB,2)],'YLim',[1 size(Hs.RGB,1)]);






%================================================
% * Callback functions for |uicontrol()| elements
%================================================





function nextCellB_Callback(hObject, eventdata)
% Button/key press to save the current set of segmentations to the
% data***.mat file and then proceed to next set of image files.
Hs = guidata(gcbo);

% NEXT/SAVE button is clicked. Disable buttons to avoid double input
set(Hs.btnHs,'Enable','off','BackgroundColor','white');
set(Hs.nextCellB,'String','Loading...');
drawnow;


load(Hs.dataFiles(Hs.currentFile).name);
objects(Hs.currentObject).metadata.type = get(Hs.typeBtns(Hs.currentSelection), 'String');
save(Hs.dataFiles(Hs.currentFile).name, 'objects')
clear objects;
Hs.currObjs = [];
load(Hs.dataFiles(Hs.currentFile).name);
if Hs.currentObject == numel(objects)
    Hs.currentObject = 1;
    Hs.currentFile = Hs.currentFile+1;
    if Hs.currentFile > numel(Hs.dataFiles)
        delete(gcf)
        return;
    else
    end
    load(Hs.dataFiles(Hs.currentFile).name)
    while isempty(objects)
        Hs.currentFile = Hs.currentFile+1;
        load(Hs.dataFiles(Hs.currentFile).name)
    end
    
else
    Hs.currentObject = Hs.currentObject+1;
end
if Hs.currentFile == numel(Hs.dataFiles) && Hs.currentObject == numel(objects)
    set(Hs.nextCellB,'String','Save & Quit');
else
    set(Hs.nextCellB,'String','Save & Next');
end


Hs = makeAndShowOverlay(Hs);
% NEXT/SAVE button actions done. Enable buttons
set(Hs.btnHs,'Enable','on','BackgroundColor','factory');
drawnow;
Hs = setFocusToFigure(Hs); % remove focus from uicontrol
guidata(gcbo, Hs);




function prevFileB_Callback(hObject, eventdata)
% Button/key press to save the current set of segmentations to the
% data***.mat file and then proceed to next set of image files.
Hs = guidata(gcbo);

% NEXT/SAVE button is clicked. Disable buttons to avoid double input
set(Hs.btnHs,'Enable','off','BackgroundColor','white');
set(Hs.nextCellB,'String','Loading...');
drawnow;

load(Hs.dataFiles(Hs.currentFile).name);
objects(Hs.currentObject).metadata.type = get(Hs.typeBtns(Hs.currentSelection), 'String');
save(Hs.dataFiles(Hs.currentFile).name, 'objects')
% Save the |image_object|s to data***.mat file
clear objects;
Hs.currObjs = [];
load(Hs.dataFiles(Hs.currentFile).name);
if Hs.currentObject == 1
    if Hs.currentFile == 1
        fprintf(1,'NOTICE: Already at the first cell\n');
        set(Hs.btnHs,'Enable','on','BackgroundColor','factory');
        drawnow;
        Hs = setFocusToFigure(Hs); % remove focus from uicontrol
        guidata(gcbo, Hs);
        return;
    end
    Hs.currentFile = Hs.currentFile-1;
    load(Hs.dataFiles(Hs.currentFile).name)
    while isempty(objects)
        Hs.currentFile = Hs.currentFile-1;
        load(Hs.dataFiles(Hs.currentFile).name)
    end
    Hs.currentObject = numel(objects);
else
    Hs.currentObject = Hs.currentObject - 1;
end

Hs = makeAndShowOverlay(Hs);
% NEXT/SAVE button actions done. Enable buttons
set(Hs.btnHs,'Enable','on','BackgroundColor','factory');
drawnow;
Hs = setFocusToFigure(Hs); % remove focus from uicontrol
guidata(gcbo, Hs);


function HeteroBtn_Callback(hObject,eventdata)
%Selection call back for heterokaryon
Hs = guidata(gcbo);
Hs.currentSelection = 1;
set(Hs.HeteroBtn, 'Value', 1);
set(Hs.GMBtn, 'Value', 0);
set(Hs.CRLBtn, 'Value', 0);
set(Hs.PolyBtn, 'Value', 0);
set(Hs.HomoBtn, 'Value', 0);
Hs = setFocusToFigure(Hs); % remove focus from uicontrol
guidata(gcbo,Hs);


function GMBtn_Callback(hObject,eventdata)
%Selection call back for GM cell
Hs = guidata(gcbo);
Hs.currentSelection = 2;
set(Hs.HeteroBtn, 'Value', 0);
set(Hs.GMBtn, 'Value', 1);
set(Hs.CRLBtn, 'Value', 0);
set(Hs.PolyBtn, 'Value', 0);
set(Hs.HomoBtn, 'Value', 0);
Hs = setFocusToFigure(Hs); % remove focus from uicontrol
guidata(gcbo,Hs);

function CRLBtn_Callback(hObject,eventdata)
%Selection call back for CRL cells
Hs = guidata(gcbo);
Hs.currentSelection = 3;
set(Hs.HeteroBtn, 'Value', 0);
set(Hs.GMBtn, 'Value', 0);
set(Hs.CRLBtn, 'Value', 1);
set(Hs.PolyBtn, 'Value', 0);
set(Hs.HomoBtn, 'Value', 0);
Hs = setFocusToFigure(Hs); % remove focus from uicontrol
guidata(gcbo,Hs);

function PolyBtn_Callback(hObject,eventdata)
%Selection call back for Polynucleated cells
Hs = guidata(gcbo);
Hs.currentSelection = 4;
set(Hs.HeteroBtn, 'Value', 0);
set(Hs.GMBtn, 'Value', 0);
set(Hs.CRLBtn, 'Value', 0);
set(Hs.PolyBtn, 'Value', 1);
set(Hs.HomoBtn, 'Value', 0);
Hs = setFocusToFigure(Hs); % remove focus from uicontrol
guidata(gcbo,Hs);

function HomoBtn_Callback(hObject,eventdata)
%Selection call back for homokaryons
Hs = guidata(gcbo);
Hs.currentSelection = 5;
set(Hs.HeteroBtn, 'Value', 0);
set(Hs.GMBtn, 'Value', 0);
set(Hs.CRLBtn, 'Value', 0);
set(Hs.PolyBtn, 'Value', 0);
set(Hs.HomoBtn, 'Value', 1);
Hs = setFocusToFigure(Hs); % remove focus from uicontrol
guidata(gcbo,Hs);



% WARNING!!!!  JAVAFRAME WILL BE DEPRECATED SOON!?!?! Who knows when but
% this is the only viable option to return focus to the main figure after
% using one of the uicontrol objects. We need focus on the figure so we can
% use the figure1_KeyReleaseFcn easily without requiring clicking on the
% figure to get keyboard shortcuts
function Hs = setFocusToFigure(Hs)
warning off MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame
javaFrame = get(Hs.figH,'JavaFrame');
javaFrame.getAxisComponent.requestFocus;


function KeyPressCallback(src,evnt)
%KeyPressFcn automatically takes in two inputs.
%Src is the object that was active when the keypress occurred.
%Evnt stores the data for the key pressed
Hs = guidata(src);
%Brings in the handles structure in to the function.

k = evnt.Key; %k is the key that is pressed.

if strcmp(k,'uparrow')
    pause(0.001)
    
    Hs.currentSelection = Hs.currentSelection-1;
    if Hs.currentSelection == 0
        Hs.currentSelection = 5;
    end
    switch Hs.currentSelection
        case 1
            hObject = Hs.HeteroBtn;
            HeteroBtn_Callback(hObject, []);
        case 2
            hObject = Hs.GMBtn;
            GMBtn_Callback(hObject, []);
        case 3
            hObject = Hs.CRLBtn;
            CRLBtn_Callback(hObject, []);
        case 4
            hObject = Hs.PolyBtn;
            PolyBtn_Callback(hObject, []);
        case 5
            hObject = Hs.HomoBtn;
            HomoBtn_Callback(hObject, []);
    end
    %Do the same thing to map to all the other callbacks.
elseif strcmp(k,'downarrow')
    pause(0.001)
    Hs.currentSelection = Hs.currentSelection+1;
    if Hs.currentSelection == 6
        Hs.currentSelection = 1;
    end
    switch Hs.currentSelection
        case 1
            hObject = Hs.HeteroBtn;
            HeteroBtn_Callback(hObject, []);
        case 2
            hObject = Hs.GMBtn;
            GMBtn_Callback(hObject, []);
        case 3
            hObject = Hs.CRLBtn;
            CRLBtn_Callback(hObject, []);
        case 4
            hObject = Hs.PolyBtn;
            PolyBtn_Callback(hObject, []);
        case 5
            hObject = Hs.HomoBtn;
            HomoBtn_Callback(hObject, []);
    end
    
elseif strcmp(k,'rightarrow')
    pause(0.001)
    hObject=Hs.nextCellB;
    nextCellB_Callback(hObject,[]);
elseif strcmp(k,'leftarrow')
    pause(0.001)
    hObject=Hs.prevCellB;
    prevFileB_Callback(hObject,[]);
end


function [Hs] = buildGUI(Hs,figH)
Hs.imgAx = axes('Parent',figH,...
    'Units','normalized',...
    'Position',[0.009,0.012,.682,.97],...
    'XTick',[],'YTick',[]);

Hs.HeteroBtn  = uicontrol('Parent',figH,...
    'Style','radiobutton',...
    'Callback',@HeteroBtn_Callback,...
    'Units','normalized',...
    'Position',[0.762 0.600 0.2 0.041],...
    'String','Heterokaryon','FontSize',12,...
    'Value',1,'ForegroundColor',[.1 .1 .1]);
Hs.GMBtn = uicontrol('Parent',figH,...
    'Style','radiobutton',...
    'Callback',@GMBtn_Callback,...
    'Units','normalized',...
    'Position',[0.762 0.550 0.2 0.041],...
    'String','GM cell','FontSize',12,...
    'Value',0,'ForegroundColor',[.1 .1 .1]);
Hs.CRLBtn  = uicontrol('Parent',figH,...
    'Style','radiobutton',...
    'Callback',@CRLBtn_Callback,...
    'Units','normalized',...
    'Position',[0.762 0.500 0.2 0.041],...
    'String','CRL cell','FontSize',12,...
    'Value',0,'ForegroundColor',[.1 .1 .1]);
Hs.PolyBtn = uicontrol('Parent',figH,...
    'Style','radiobutton',...
    'Callback',@PolyBtn_Callback,...
    'Units','normalized',...
    'Position',[0.762 0.450 0.2 0.041],...
    'String','Polynucleated','FontSize',12,...
    'Value',0,'ForegroundColor',[.1 .1 .1]);
Hs.HomoBtn = uicontrol('Parent',figH,...
    'Style','radiobutton',...
    'Callback',@HomoBtn_Callback,...
    'Units','normalized',...
    'Position',[0.762 0.400 0.2 0.041],...
    'String','Homokaryon','FontSize',12,...
    'Value',0,'ForegroundColor',[.1 .1 .1]);
Hs.prevCellB = uicontrol('Parent',figH,...
    'Style','pushbutton',...
    'Callback',@prevCellB_Callback,...
    'Units','normalized',...
    'Position',[0.743 0.233 0.1 0.048],...
    'String','Prev Cell');
Hs.nextCellB = uicontrol('Parent',figH,...
    'Style','pushbutton',...
    'Callback',@nextCellB_Callback,...
    'Units','normalized',...
    'Position',[0.851 0.233 0.1 0.048]);
if numel(Hs.dataFiles) == 1
    set(Hs.nextCellB,'String','Save & Quit');
    set(Hs.prevCellB,'Enable','off');
else
    set(Hs.nextCellB,'String','Save & Next');
end
set([Hs.prevCellB Hs.nextCellB],'FontSize',12,'FontWeight','bold');
Hs.btnHs = [Hs.prevCellB Hs.nextCellB];
Hs.typeBtns = [Hs.HeteroBtn, Hs.GMBtn, Hs.CRLBtn, Hs.PolyBtn, Hs.HomoBtn];
Hs.currentSelection = 1;



