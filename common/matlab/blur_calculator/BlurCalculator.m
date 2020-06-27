function varargout = BlurCalculator(varargin)
% BLURCALCULATOR MATLAB code for BlurCalculator.fig
%      BLURCALCULATOR, by itself, creates a new BLURCALCULATOR or raises the existing
%      singleton*.
%
%      H = BLURCALCULATOR returns the handle to a new BLURCALCULATOR or the handle to
%      the existing singleton*.
%
%      BLURCALCULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BLURCALCULATOR.M with the given input arguments.
%
%      BLURCALCULATOR('Property','Value',...) creates a new BLURCALCULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BlurCalculator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BlurCalculator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BlurCalculator

% Last Modified by GUIDE v2.5 19-Nov-2018 11:28:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BlurCalculator_OpeningFcn, ...
                   'gui_OutputFcn',  @BlurCalculator_OutputFcn, ...
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


% --- Executes just before BlurCalculator is made visible.
function BlurCalculator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BlurCalculator (see VARARGIN)

% Choose default command line output for BlurCalculator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BlurCalculator wait for user response (see UIRESUME)
% uiwait(handles.figure1);
updatePlot(handles)
% xlabel(handles.plot_axes,'Distance From Camera (m)', 'fontweight','bold');
% ylabel(handles.plot_axes,'Blur Radius (Pixels)', 'fontweight','bold');
% title(handles.plot_axes,'Object Distance vs. Radius of Blur', 'fontweight','bold');


% --- Outputs from this function are returned to the command line.
function varargout = BlurCalculator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function pixel_size_TB_Callback(hObject, eventdata, handles)
% hObject    handle to pixel_size_TB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_size_TB as text
%        str2double(get(hObject,'String')) returns contents of pixel_size_TB as a double
updatePlot(handles);


% --- Executes during object creation, after setting all properties.
function pixel_size_TB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_size_TB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f_num_TB_Callback(hObject, eventdata, handles)
% hObject    handle to f_num_TB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_num_TB as text
%        str2double(get(hObject,'String')) returns contents of f_num_TB as a double
updatePlot(handles);


% --- Executes during object creation, after setting all properties.
function f_num_TB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_num_TB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f_length_TB_Callback(hObject, eventdata, handles)
% hObject    handle to f_length_TB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f_length_TB as text
%        str2double(get(hObject,'String')) returns contents of f_length_TB as a double
updatePlot(handles);


% --- Executes during object creation, after setting all properties.
function f_length_TB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_length_TB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function fd_SLD_Callback(hObject, eventdata, handles)
% hObject    handle to fd_SLD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_val = floor(get(hObject,'Value')*10)/10;
set(handles.fd_TB,'String',num2str(slider_val));
updatePlot(handles);




% --- Executes during object creation, after setting all properties.
function fd_SLD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fd_SLD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function fd_TB_Callback(hObject, eventdata, handles)
% hObject    handle to fd_TB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fd_TB as text
%        str2double(get(hObject,'String')) returns contents of fd_TB as a double
fd_val = str2double(get(hObject,'String'));
set(handles.fd_SLD,'Value',fd_val);
updatePlot(handles);



% --- Executes during object creation, after setting all properties.
function fd_TB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fd_TB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function plot_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate plot_axes



function xmin_lim_Callback(hObject, eventdata, handles)
% hObject    handle to xmin_lim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmin_lim as text
%        str2double(get(hObject,'String')) returns contents of xmin_lim as a double
% x_min = str2double(get(hObject,'String'));
% curr_lim = get(handles.plot_axes,'xlim');
% curr_lim(1) = x_min;
% set(handles.plot_axes,'xlim', curr_lim);
updatePlot(handles)

% --- Executes during object creation, after setting all properties.
function xmin_lim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmin_lim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xmax_lim_Callback(hObject, eventdata, handles)
% hObject    handle to xmax_lim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmax_lim as text
%        str2double(get(hObject,'String')) returns contents of xmax_lim as a double
% x_max = str2double(get(hObject,'String'));
% curr_lim = get(handles.plot_axes,'xlim');
% curr_lim(2) = x_max;
% set(handles.plot_axes,'xlim', curr_lim);
updatePlot(handles);
%fprintf(1,'%f, %f\n',get(handles.plot_axes,'xlim'))

% --- Executes during object creation, after setting all properties.
function xmax_lim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmax_lim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function updatePlot(handles)

limits = [str2double(get(handles.xmin_lim,'String')),str2double(get(handles.xmax_lim,'String'))];
px = str2double(get(handles.pixel_size_TB,'String'));
f_num = str2double(get(handles.f_num_TB,'String'));
fl = str2double(get(handles.f_length_TB,'String'));
d_o = str2double(get(handles.fd_TB,'String'));

d_o = d_o * 1000;

c_lim = 1*px;

H = ((fl*fl)/(c_lim*f_num));% + fl;
% Dn = (df*1000*(fl*fl)/(fl*fl + f_num*c_lim*(df*1000-fl)))/1000;
% Df = (df*1000*(fl*fl)/(fl*fl - f_num*c_lim*(df*1000-fl)))/1000;

%Dn = ((d_o*1000)*(H-fl)/(H+(d_o*1000)-2*fl))/1000;
%Df = ((d_o*1000)*(H-fl)/(H-(d_o*1000)))/1000;
Dn = ((d_o)*(fl*fl)/(fl*fl+c_lim*f_num*(d_o-fl)))/1000;
Df = ((d_o)*(fl*fl)/(fl*fl-c_lim*f_num*(d_o-fl)))/1000;

DOF = Df-Dn;

set(handles.d_near_tb, 'String', num2str(Dn));
if(Df >= 0)
    set(handles.d_far_tb, 'String', num2str(Df));
    set(handles.dof_tb, 'String', num2str(DOF));
else
    set(handles.d_far_tb, 'String', 'Inf');
    set(handles.dof_tb, 'String', 'Inf'); 
    Df = inf;
end

% df = focal distance (m)
% fl = focal length (mm)
% f_num = lens f-number
% px = pixel size (mm)

[S_range, CoC, CoC_max] = blurCalc(f_num, fl, d_o, limits*1000);
%disp(CoC_max);

y_axis_ticks = [0:px:(CoC_max+px)];
y_axis_labels = num2str(y_axis_ticks'/px);
hold off
plot(handles.plot_axes, S_range/1000, CoC,'.-b');
hold on
plot(handles.plot_axes, S_range/1000, ceil(CoC/px)*px,'-k');
plot(handles.plot_axes, [0, S_range(end)/1000], [CoC_max, CoC_max],'-g');
stem([Dn, Df],[CoC_max, CoC_max],'.r');

set(handles.plot_axes,'Xlim', limits);
set(handles.plot_axes,'ylim', [0, CoC_max+px],'fontweight','bold')
set(handles.plot_axes, 'yTickLabel', y_axis_labels, 'YTick',y_axis_ticks);

xlabel(handles.plot_axes,'Distance From Camera (m)', 'fontweight','bold');
ylabel(handles.plot_axes,'Blur Radius (Pixels)', 'fontweight','bold');
title(handles.plot_axes,'Object Distance vs. Radius of Blur', 'fontweight','bold');
grid on



function d_near_tb_Callback(hObject, eventdata, handles)
% hObject    handle to d_near_tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_near_tb as text
%        str2double(get(hObject,'String')) returns contents of d_near_tb as a double


% --- Executes during object creation, after setting all properties.
function d_near_tb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_near_tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d_far_tb_Callback(hObject, eventdata, handles)
% hObject    handle to d_far_tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_far_tb as text
%        str2double(get(hObject,'String')) returns contents of d_far_tb as a double


% --- Executes during object creation, after setting all properties.
function d_far_tb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_far_tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dof_tb_Callback(hObject, eventdata, handles)
% hObject    handle to dof_tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dof_tb as text
%        str2double(get(hObject,'String')) returns contents of dof_tb as a double


% --- Executes during object creation, after setting all properties.
function dof_tb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dof_tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
