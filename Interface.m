function varargout = Interface(varargin)
% INTERFACE MATLAB code for Interface.fig
%      INTERFACE, by itself, creates a new INTERFACE or raises the existing
%      singleton*.
%
%      H = INTERFACE returns the handle to a new INTERFACE or the handle to
%      the existing singleton*.
%
%      INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERFACE.M with the given input arguments.
%
%      INTERFACE('Property','Value',...) creates a new INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Interface

% Last Modified by GUIDE v2.5 17-Aug-2020 17:33:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Interface_OpeningFcn, ...
                   'gui_OutputFcn',  @Interface_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Interface is made visible.
function Interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Interface (see VARARGIN)

% Choose default command line output for Interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in injury.
function injury_Callback(hObject, eventdata, handles)
% hObject    handle to injury (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns injury contents as cell array
%        contents{get(hObject,'Value')} returns selected item from injury
ListBoxString = get(handles.injury,'String');
ListBoxValue = get(handles.injury,'Value');
global folderName
if ListBoxValue==1 
    folderName = 'Recruitment/huttenlocher_injury/';
         cc = 0.99; % scaling coefficient for the scale bar
         T  = 0.5;  % time increment
         ListBoxContents = {1,2,3};
         nFish = 3; % number of fish to process
elseif ListBoxValue==2
    folderName = 'Recruitment/mild_injury/';
         cc = 0.93; % scaling coefficient for the scale bar
         T  = 1.5;  % time increment
         ListBoxContents = {1,2,3,4,5};
         nFish = 5; % number of fish to process
elseif ListBoxValue==3
    folderName = 'Recruitment/normal_injury/';
         cc = 1.71; % scaling coefficient for the scale bar
         T  = 2;    % time increment
         ListBoxContents = {1,2,3,4,5,6};
         nFish = 6; % number of fish to process
elseif ListBoxValue==4
    folderName = 'Recruitment/severe_injury/';
         cc = 1; % scaling coefficient for the scale bar
         T  = 2; % time increment
         ListBoxContents = {1,2};
         nFish = 2; % number of fish to process
end 
set(handles.data, 'String',ListBoxContents);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function injury_CreateFcn(hObject, eventdata, handles)
% hObject    handle to injury (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in data.
function data_Callback(hObject, eventdata, handles)
% hObject    handle to data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ListBoxString = get(handles.data,'String');
ListBoxValue = get(handles.data,'Value');
global iFish
if ListBoxValue==1
    iFish = 1;
elseif ListBoxValue==2
    iFish = 2;
elseif ListBoxValue==3
    iFish = 3;
elseif ListBoxValue==4
    iFish = 4;
elseif ListBoxValue==5
    iFish = 5;
elseif ListBoxValue==6
    iFish = 6;
end
% Hints: contents = cellstr(get(hObject,'String')) returns data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data


% --- Executes during object creation, after setting all properties.
function data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in figure.
function figure_Callback(hObject, eventdata, handles, fold)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ListBoxString = get(handles.figure,'String');
ListBoxValue = get(handles.figure,'Value');
global folderName
global iFish
if ListBoxValue==1
    % Load tracking data
    load([folderName,'tracks_' num2str(iFish)]);
    % Load brightfield image
    A = imread([folderName,'bf_' num2str(iFish),'.png']); % for huttenlocher injury 
    Cnt = rgb2gray(A); 
    [y_max,x_max,zz] = size(A);
    BW = fishmask(Cnt);
    AA = double(BW); % create a surface
    Xx = 1:1:size(A,2); % create the grid of x coords
    Yy = 1:1:size(A,1); % create the grid of y coords
    [Yy_grid,Xx_grid] = meshgrid(Xx,Yy); % mesh
    white=[1,1,1]; % surface colour
    gg = [0.8,0.8,0.8]; % extra colour for cells
    set(gcf,'color','w');
    axes(handles.axes1); % Switch current axes to axes1
    imshow(A);
    
elseif ListBoxValue == 2
    load([folderName,'tracks_' num2str(iFish)]);
% Load brightfield image
A = imread([folderName,'bf_' num2str(iFish),'.png']); % for huttenlocher injury 
Cnt = rgb2gray(A); 
[y_max,x_max,zz] = size(A);
% Load the mask
BW = fishmask(Cnt);
% For data split into hours
hour = 0;
if hour~=0
    X = XX;
    Y = YY;
    nTracks = nTracks_t;
end
Tracks = 1:nTracks;
% Create a frame around the image to extend basis function support;
padH = 0; % vertical padding
padW = 200; % horizontal padding
A = padarray(A,[padH, padW]); % creating a padded image
for i=Tracks
    X{i}(:,1) = X{i}(:,1) + padW;
    X{i}(:,2) = X{i}(:,2) + padH;
    Y{i}(:,1) = Y{i}(:,1) + padW;
    Y{i}(:,2) = Y{i}(:,2) + padH;
end
[y_max,x_max,zz] = size(A);
x_lim = [padW x_max-padW];
y_lim = [padH y_max-padH];

   AA = double(BW); % create a surface
Xx = 1:1:size(A,2); % create the grid of x coords
Yy = 1:1:size(A,1); % create the grid of y coords
[Yy_grid,Xx_grid] = meshgrid(Xx,Yy); % mesh
white=[1,1,1]; % surface colour
gg = [0.8,0.8,0.8]; % extra colour for cells 
    
   
axes(handles.axes1);
imshow(A); hold on;
% hold on;
for j = Tracks
   plot(Y{j}(:,1),Y{j}(:,2),'-w','LineWidth',1); hold on; 
end
surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
view(2)
xlim(x_lim);ylim(y_lim);
hold on;
line([250,250+100*cc],[y_max-20,y_max-20],[2,2],'Color','k','LineWidth',5);
txt = ('100 $\mu$m');
text(250,y_max-70, 2,txt,'Color','k','FontSize',20)
set(gca,'Ydir','reverse')

end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
