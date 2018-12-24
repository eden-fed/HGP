function varargout = LoadMesh(varargin)
% LOADMESH MATLAB code for LoadMesh.fig
%      LOADMESH, by itself, creates a new LOADMESH or raises the existing
%      singleton*.
%
%      H = LOADMESH returns the handle to a new LOADMESH or the handle to
%      the existing singleton*.
%
%      LOADMESH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADMESH.M with the given input arguments.
%
%      LOADMESH('Property','Value',...) creates a new LOADMESH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LoadMesh_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LoadMesh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LoadMesh

% Last Modified by GUIDE v2.5 29-Nov-2018 17:08:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LoadMesh_OpeningFcn, ...
                   'gui_OutputFcn',  @LoadMesh_OutputFcn, ...
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


% --- Executes just before LoadMesh is made visible.
function LoadMesh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LoadMesh (see VARARGIN)

% Choose default command line output for LoadMesh
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LoadMesh wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LoadMesh_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function obj_Callback(hObject, eventdata, handles)
% hObject    handle to obj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of obj as text
%        str2double(get(hObject,'String')) returns contents of obj as a double


% --- Executes during object creation, after setting all properties.
function obj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to obj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vf_Callback(hObject, eventdata, handles)
% hObject    handle to vf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vf as text
%        str2double(get(hObject,'String')) returns contents of vf as a double


% --- Executes during object creation, after setting all properties.
function vf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in objSelect.
function objSelect_Callback(hObject, eventdata, handles)
% hObject    handle to objSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FilterSpec = '*.obj';
[FileName,PathName] = uigetfile(FilterSpec);
str = [PathName,FileName];
set(handles.obj,'string',str);
str(end)=[];
str(end)=[];
str(end)=[];
str=[str,'ffield'];
set(handles.vf,'string',str);

% --- Executes on button press in vfSelect.
function vfSelect_Callback(hObject, eventdata, handles)
% hObject    handle to vfSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FilterSpec = {'*.ffield';'*.mat'};
[FileName,PathName] = uigetfile(FilterSpec);
str = [PathName,FileName];
set(handles.vf,'string',str);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
objLocation = get(handles.obj,'string');
objLocation = regexprep(objLocation,'\','\\\');
vfLocation = get(handles.vf,'string');
vfLocation = regexprep(vfLocation,'\','\\\');
methodIndex = get(handles.method, 'value');

assignin('base','objLocation',objLocation);
assignin('base','vfLocation',vfLocation);
assignin('base','methodIndex',methodIndex);

delete(gcbf)


% --- Executes on selection change in method.
function method_Callback(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method




% --- Executes during object creation, after setting all properties.
function method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
