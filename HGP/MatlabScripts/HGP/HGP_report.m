function varargout = HGP_report(varargin)
% HGP_REPORT MATLAB code for HGP_report.fig
%      HGP_REPORT, by itself, creates a new HGP_REPORT or raises the existing
%      singleton*.
%
%      H = HGP_REPORT returns the handle to a new HGP_REPORT or the handle to
%      the existing singleton*.
%
%      HGP_REPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HGP_REPORT.M with the given input arguments.
%
%      HGP_REPORT('Property','Value',...) creates a new HGP_REPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HGP_report_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HGP_report_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HGP_report

% Last Modified by GUIDE v2.5 20-Nov-2017 13:24:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HGP_report_OpeningFcn, ...
                   'gui_OutputFcn',  @HGP_report_OutputFcn, ...
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


% --- Executes just before HGP_report is made visible.
function HGP_report_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HGP_report (see VARARGIN)

% Choose default command line output for HGP_report
handles.output = hObject;

% Update handles structure
HGP = evalin('base','HGP');
visualizeFlag = evalin('base','visMatlab');
if visualizeFlag == 0
   set(handles.visualize,'Enable','off'); 
end
set(handles.numV, 'String', num2str(length(HGP.V)));
set(handles.numF, 'String', num2str(length(HGP.F)));
set(handles.numC, 'String', num2str(length(HGP.conesIndices)));
set(handles.numB, 'String', num2str(HGP.Result.numBorders));
set(handles.numG, 'String', num2str(HGP.Result.genus));

set(handles.itNum, 'String', num2str(length(HGP.Result.timeVector)));
%useFix
if HGP.FrameFix.use == 1
    set(handles.useFix, 'String', 'Yes');
end

set(handles.kMin, 'String', num2str(fix(min(HGP.Result.k)*10000000)/10000000));
set(handles.kMax, 'String', num2str(fix(max(HGP.Result.k)*10000000)/10000000));
set(handles.kAvg, 'String', num2str( fix(sum(HGP.Result.k)/length(HGP.Result.k)*10000000)/10000000) );

set(handles.itTimes, 'Data', [(1:length(HGP.Result.timeVector))' HGP.Result.timeVector'])
set(handles.itTimes, 'ColumnName',{'Iteration'; 'Time'})
set(handles.itTimes, 'ColumnWidth',{50 'auto'})

set(handles.foldOversTable, 'Data', HGP.Result.flips)
set(handles.foldOversTable, 'ColumnName',{'Triangle Index'})

set(handles.coneTable, 'Data', [HGP.Result.problemVerticesIndices HGP.Result.problemVerticesAngles HGP.Result.desiredVerticesAngles])
set(handles.coneTable, 'ColumnName',{'Index'; 'Angle'; 'Desired Angle'})
set(handles.coneTable, 'ColumnWidth',{50 'auto'})

handles.HGP = HGP;
guidata(hObject, handles);

% UIWAIT makes HGP_report wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HGP_report_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in saveLog.
function saveLog_Callback(hObject, eventdata, handles)
% hObject    handle to saveLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[name,path] = uiputfile([handles.HGP.meshName '_HGP_log.txt']);
if name == 0
    return
end
fileID = fopen([path name],'w+t');

fprintf(fileID,'*******HGP log******\n');
fprintf(fileID,'Model name: %s\n\n',handles.HGP.meshName);
fprintf(fileID,'Number of vertices: %i\n',length(handles.HGP.V));
fprintf(fileID,'Number of faces: %i\n',length(handles.HGP.F));
fprintf(fileID,'Number of cones: %i\n',length(handles.HGP.conesIndices));
fprintf(fileID,'Number of borders: %i\n',handles.HGP.Result.numBorders);
fprintf(fileID,'Genus: %i\n\n',handles.HGP.Result.genus);

fprintf(fileID,'Total number of iterations: %i\n',length(handles.HGP.Result.timeVector));
for i=1:length(handles.HGP.Result.timeVector)
    fprintf(fileID,'Iteration number %i time: %5.5f seconds\n',i,handles.HGP.Result.timeVector(i));
end

if handles.HGP.FrameFix.use == 0
    fprintf(fileID,'\nUse frame fix: No\n');
else
    fprintf(fileID,'\nUse frame fix: Yes\n');
end
    
if isempty(handles.HGP.Result.flips) == 1
    fprintf(fileID,'\nNo fold-overs!\n');
else
    fprintf(fileID,'\n**ERROR**There are %i fold-overs. Triangles indices:\n',length(handles.HGP.Result.flips));
    for i = 1:length(handles.HGP.Result.flips)
        fprintf(fileID,'%i,',handles.HGP.Result.flips(i));
    end
    fprintf(fileID,'\n');
end

if isempty(handles.HGP.Result.problemVerticesIndices) == 1
    fprintf(fileID,'All angles are o.k!\n');
else
    fprintf(fileID,'**WARNING** There are problems with the angles!:\n');
end

fprintf(fileID,'\nk min: %5.7f\n',min(handles.HGP.Result.k));
fprintf(fileID,'k max: %5.7f\n',max(handles.HGP.Result.k));
fprintf(fileID,'k avg: %5.7f\n',sum(handles.HGP.Result.k)/length(handles.HGP.Result.k));

fclose(fileID);

% --- Executes on button press in visualize.
function visualize_Callback(hObject, eventdata, handles)
% hObject    handle to visualize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

choice = questdlg('This operation might take some time. Do you want to continue?', 'Warning', ...
	'Yes','No','Yes');
if strcmp(choice,'No') == 1
    return;
end

visualizeMesh( handles.HGP );

% --- Executes on button press in exportButton.
function exportButton_Callback(hObject, eventdata, handles)
% hObject    handle to exportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name,path] = uiputfile([handles.HGP.meshName '.obj']);
if name == 0
    return
end

saveObj( path, name, handles.HGP )