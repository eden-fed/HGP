function varargout = FastHGP_settings(varargin)
% FASTHGP_SETTINGS MATLAB code for FastHGP_settings.fig
%      FASTHGP_SETTINGS, by itself, creates a new FASTHGP_SETTINGS or raises the existing
%      singleton*.
%
%      H = FASTHGP_SETTINGS returns the handle to a new FASTHGP_SETTINGS or the handle to
%      the existing singleton*.
%
%      FASTHGP_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FASTHGP_SETTINGS.M with the given input arguments.
%
%      FASTHGP_SETTINGS('Property','Value',...) creates a new FASTHGP_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FastHGP_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FastHGP_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FastHGP_settings

% Last Modified by GUIDE v2.5 23-Dec-2018 14:50:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FastHGP_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @FastHGP_settings_OutputFcn, ...
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


% --- Executes just before FastHGP_settings is made visible.
function FastHGP_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FastHGP_settings (see VARARGIN)

% Choose default command line output for FastHGP_settings
handles.output = hObject;

% Update handles structure
numV = evalin('base','numV');
numF = evalin('base','numF');
numC = evalin('base','numC');
numB = evalin('base','FastHGP.numB');

set(handles.numV, 'String', num2str(numV));
set(handles.numF, 'String', num2str(numF));
set(handles.numC, 'String', num2str(numC));
set(handles.numB, 'String', num2str(numB));

guidata(hObject, handles);

% UIWAIT makes FastHGP_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FastHGP_settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;



function segSize_Callback(hObject, eventdata, handles)
% hObject    handle to segSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of segSize as text
%        str2double(get(hObject,'String')) returns contents of segSize as a double


% --- Executes during object creation, after setting all properties.
function segSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fixCot.
function fixCot_Callback(hObject, eventdata, handles)
% hObject    handle to fixCot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fixCot


% --- Executes on button press in continueButton.
function continueButton_Callback(hObject, eventdata, handles)
% hObject    handle to continueButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
visMatlab = get(handles.visMatlab,'value');
fixCot = get(handles.fixCot,'value');
segSize = str2double(get(handles.segSize,'String'));

assignin('base','visMatlab',visMatlab);
assignin('base','fixCot',fixCot);
assignin('base','segSize',segSize);

delete(gcbf)

% --- Executes on button press in visMatlab.
function visMatlab_Callback(hObject, eventdata, handles)
% hObject    handle to visMatlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of visMatlab
