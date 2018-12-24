function varargout = FastHGP_report(varargin)
% FASTHGP_REPORT MATLAB code for FastHGP_report.fig
%      FASTHGP_REPORT, by itself, creates a new FASTHGP_REPORT or raises the existing
%      singleton*.
%
%      H = FASTHGP_REPORT returns the handle to a new FASTHGP_REPORT or the handle to
%      the existing singleton*.
%
%      FASTHGP_REPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FASTHGP_REPORT.M with the given input arguments.
%
%      FASTHGP_REPORT('Property','Value',...) creates a new FASTHGP_REPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FastHGP_report_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FastHGP_report_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FastHGP_report

% Last Modified by GUIDE v2.5 23-Dec-2018 14:49:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FastHGP_report_OpeningFcn, ...
                   'gui_OutputFcn',  @FastHGP_report_OutputFcn, ...
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


% --- Executes just before FastHGP_report is made visible.
function FastHGP_report_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FastHGP_report (see VARARGIN)

% Choose default command line output for FastHGP_report
handles.output = hObject;

% Update handles structure
FastHGP = evalin('base','FastHGP');
visualizeFlag = evalin('base','visMatlab');
if visualizeFlag == 0
   set(handles.visualize,'Enable','off'); 
end
set(handles.numV, 'String', num2str(length(FastHGP.V)));
set(handles.numF, 'String', num2str(length(FastHGP.F)));
set(handles.numC, 'String', num2str(length(FastHGP.conesIndices)));
set(handles.numB, 'String', num2str(FastHGP.numB));
set(handles.numDOF, 'String', num2str(FastHGP.numDOF));
set(handles.totalTime, 'String', num2str(FastHGP.Result.totalTime));

set(handles.kMin, 'String', num2str(fix(min(FastHGP.Result.k)*10000000)/10000000));
set(handles.kMax, 'String', num2str(fix(max(FastHGP.Result.k)*10000000)/10000000));
set(handles.kAvg, 'String', num2str( fix(sum(FastHGP.Result.k)/length(FastHGP.Result.k)*10000000)/10000000) );

set(handles.foldOversTable, 'Data', FastHGP.Result.flips)
set(handles.foldOversTable, 'ColumnName',{'Triangle Index'})

set(handles.coneTable, 'Data', [FastHGP.Result.problemVerticesIndices FastHGP.Result.problemVerticesAngles FastHGP.Result.desiredVerticesAngles])
set(handles.coneTable, 'ColumnName',{'Index'; 'Angle'; 'Desired Angle'})
set(handles.coneTable, 'ColumnWidth',{50 'auto'})

handles.FastHGP = FastHGP;
guidata(hObject, handles);

% UIWAIT makes FastHGP_report wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FastHGP_report_OutputFcn(hObject, eventdata, handles) 
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

[name,path] = uiputfile([handles.FastHGP.meshName '_FastHGP_log.txt']);
if name == 0
    return
end
fileID = fopen([path name],'w+t');

fprintf(fileID,'*******FastHGP log******\n');
fprintf(fileID,'Model name: %s\n\n',handles.FastHGP.meshName);
fprintf(fileID,'Number of vertices: %i\n',length(handles.FastHGP.V));
fprintf(fileID,'Number of faces: %i\n',length(handles.FastHGP.F));
fprintf(fileID,'Number of cones: %i\n',length(handles.FastHGP.conesIndices));
fprintf(fileID,'Number of border vertices: %i\n',handles.FastHGP.numB);
fprintf(fileID,'Number of DOF: %i\n',handles.FastHGP.numDOF);

fprintf(fileID,'\nTotal time of algorithm: %i\n',handles.FastHGP.Result.totalTime);

if isempty(handles.FastHGP.Result.flips) == 1
    fprintf(fileID,'\nNo fold-overs!\n');
else
    fprintf(fileID,'\n**ERROR**There are %i fold-overs. Triangles indices:\n',length(handles.FastHGP.Result.flips));
    for i = 1:length(handles.FastHGP.Result.flips)
        fprintf(fileID,'%i,',handles.FastHGP.Result.flips(i));
    end
    fprintf(fileID,'\n');
end

if isempty(handles.FastHGP.Result.problemVerticesIndices) == 1
    fprintf(fileID,'All angles are o.k!\n');
else
    fprintf(fileID,'**WARNING** There are problems with the angles!:\n');
end

fprintf(fileID,'\nk min: %5.7f\n',min(handles.FastHGP.Result.k));
fprintf(fileID,'k max: %5.7f\n',max(handles.FastHGP.Result.k));
fprintf(fileID,'k avg: %5.7f\n',sum(handles.FastHGP.Result.k)/length(handles.FastHGP.Result.k));

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

visualizeMesh( handles.FastHGP );

% --- Executes on button press in exportButton.
function exportButton_Callback(hObject, eventdata, handles)
% hObject    handle to exportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name,path] = uiputfile([handles.FastHGP.meshName '.obj']);
if name == 0
    return
end

saveObj( path, name, handles.FastHGP )
