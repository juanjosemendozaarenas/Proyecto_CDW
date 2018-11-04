function varargout = Interfaz(varargin)
% INTERFAZ MATLAB code for Interfaz.fig
%      INTERFAZ, by itself, creates a new INTERFAZ or raises the existing
%      singleton*.
%
%      H = INTERFAZ returns the handle to a new INTERFAZ or the handle to
%      the existing singleton*.
%
%      INTERFAZ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERFAZ.M with the given input arguments.
%
%      INTERFAZ('Property','Value',...) creates a new INTERFAZ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Interfaz_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Interfaz_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Interfaz

% Last Modified by GUIDE v2.5 03-Aug-2018 10:49:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Interfaz_OpeningFcn, ...
                   'gui_OutputFcn',  @Interfaz_OutputFcn, ...
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


% --- Executes just before Interfaz is made visible.
function Interfaz_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Interfaz (see VARARGIN)

clc;

% Choose default command line output for Interfaz
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Interfaz wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Interfaz_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Jobs.
function Jobs_Callback(hObject, eventdata, handles)
% hObject    handle to Jobs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Jobs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Jobs

str = get(hObject,'String');
val = get(hObject,'Value');

[handles.DBL, handles.SZ, handles.NUP, handles.T_DBL,...
            handles.T_SZ, handles.T_NUP, handles.T_NDN ,...
            handles.MPn, handles.MnP, handles.MPt,...
            handles.MtP, handles.param] = get_data(str2double(str{val}));

switch str{val}
    case '1'
        handles.F = 'U = -2';
        handles.P = 'V';
    case '2'
        handles.F = 'U = 1';
        handles.P = 'V';
    case '3'
        handles.F = 'V = 4';
        handles.P = 'U';
    case '96'
        handles.F = 'V = 4';
        handles.P = 'U';
    case '4'
        handles.F = 'V = -0.2';
        handles.P = 'U';
end

% Save handles data
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Jobs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Jobs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exp_dbl.
function exp_dbl_Callback(hObject, eventdata, handles)
% hObject    handle to exp_dbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(1)
surf(handles.MPn,handles.MnP,handles.DBL);
title(['Doble ocupación para ' handles.F])
ylabel('Sitio')
xlabel(handles.P)
zlabel('Doble Ocupación')

% --- Executes on button press in exp_nup.
function exp_nup_Callback(hObject, eventdata, handles)
% hObject    handle to exp_nup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(2)
surf(handles.MPn,handles.MnP,handles.NUP);
title(['Espines arriba para ' handles.F])
ylabel('Sitio')
xlabel(handles.P)
zlabel('Espines arriba')

% --- Executes on button press in exp_sz.
function exp_sz_Callback(hObject, eventdata, handles)
% hObject    handle to exp_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(3)
surf(handles.MPn,handles.MnP,handles.SZ);
title(['Magnetización para ' handles.F])
ylabel('Sitio')
xlabel(handles.P)
zlabel('Magnetización')

% --- Executes on button press in tc_dbl.
function tc_dbl_Callback(hObject, eventdata, handles)
% hObject    handle to tc_dbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(4)
surf(handles.MPt,handles.MtP,handles.T_DBL);
title(['Correlaciones temporales DBL para ' handles.F])
ylabel('Tiempo')
xlabel(handles.P)
zlabel('TC')

% figure()
% ax1 = subplot(3,1,1);
% plot(handles.V,handles.plano(chosen_times(1),:),'b','linewidth',2)
% title(['Time Correlations at times ' num2str(handles.T(chosen_times(1)))...
%     ' and ' num2str(handles.T(chosen_times(2))) ' for '  handles.sim])
% hold on
% plot(handles.V,handles.plano(chosen_times(2),:),'r','linewidth',2)
% legend(['T=' num2str(handles.T(chosen_times(1)))],...
%     ['T=' num2str(handles.T(chosen_times(2)))],'Location', 'Best')
% hold off
% 
% ax2 = subplot(3,1,2);
% plot(handles.V,D(chosen_times(1),:),'g','linewidth',2)
% title('First Derivative')
% 
% ax3 = subplot(3,1,3);
% plot(handles.V,D2(chosen_times(1),:),'m','linewidth',2)
% title('Second Derivative')
% xlabel(handles.xaxis)
% linkaxes([ax1,ax2,ax3],'x')


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Cerrar.
function Cerrar_Callback(hObject, eventdata, handles)
% hObject    handle to Cerrar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all



function tiempo_Callback(hObject, eventdata, handles)
% hObject    handle to tiempo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tiempo as text
%        str2double(get(hObject,'String')) returns contents of tiempo as a double


% --- Executes during object creation, after setting all properties.
function tiempo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tiempo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tc_sz.
function tc_sz_Callback(hObject, eventdata, handles)
% hObject    handle to tc_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(5)
surf(handles.MPt,handles.MtP,handles.T_SZ);
title(['Correlaciones temporales SZ para ' handles.F])
ylabel('Tiempo')
xlabel(handles.P)
zlabel('TC')


% --- Executes on button press in tc_nup.
function tc_nup_Callback(hObject, eventdata, handles)
% hObject    handle to tc_nup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(6)
surf(handles.MPt,handles.MtP,handles.T_NUP);
title(['Correlaciones temporales NUP para ' handles.F])
ylabel('Tiempo')
xlabel(handles.P)
zlabel('TC')


% --- Executes on button press in tc_n.
function tc_n_Callback(hObject, eventdata, handles)
% hObject    handle to tc_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(7)
surf(handles.MPt,handles.MtP,0.5*(handles.T_NUP+handles.T_NDN));
title(['Correlaciones temporales N para ' handles.F])
ylabel('Tiempo')
xlabel(handles.P)
zlabel('TC')


% --- Executes on button press in lg_dbl.
function lg_dbl_Callback(hObject, eventdata, handles)
% hObject    handle to lg_dbl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LG, MtP, MPt] = Leggett_Garg(handles.T_DBL, handles.param);
figure(8)
surf(MtP,MPt,LG);
title(['Desigualdades de Leggett-Garg en DBL para ' handles.F])
ylabel(handles.P)
xlabel('\tau')
zlabel('LGI')


% --- Executes on button press in lg_sz.
function lg_sz_Callback(hObject, eventdata, handles)
% hObject    handle to lg_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LG, MtP, MPt] = Leggett_Garg(handles.T_SZ, handles.param);
figure(9)
surf(MtP,MPt,LG);
title(['Desigualdades de Leggett-Garg en SZ para ' handles.F])
ylabel(handles.P)
xlabel('\tau')
zlabel('LGI')


% --- Executes on button press in lg_nup.
function lg_nup_Callback(hObject, eventdata, handles)
% hObject    handle to lg_nup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LG, MtP, MPt] = Leggett_Garg(handles.T_NUP, handles.param);
figure(10)
surf(MtP,MPt,LG);
title(['Desigualdades de Leggett-Garg en NUP para ' handles.F])
ylabel(handles.P)
xlabel('\tau')
zlabel('LGI')


% --- Executes on button press in lg_n.
function lg_n_Callback(hObject, eventdata, handles)
% hObject    handle to lg_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[LG, MtP, MPt] = Leggett_Garg(0.5*(handles.T_NUP+handles.T_NDN), handles.param);
figure(11)
surf(MtP,MPt,LG);
title(['Desigualdades de Leggett-Garg en N para ' handles.F])
ylabel(handles.P)
xlabel('\tau')
zlabel('LGI')
