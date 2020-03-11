function varargout = INSPR_ast_GUI(varargin)
% INSPR_ast_GUI MATLAB code for INSPR_ast_GUI.fig
% This code is used for astigmatism based setup
%
% (C) Copyright 2019                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, December 2019
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @INSPR_ast_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @INSPR_ast_GUI_OutputFcn, ...
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





% --- Executes just before INSPR_ast_GUI is made visible.
function INSPR_ast_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to INSPR_ast_GUI (see VARARGIN)

% Choose default command line output for INSPR_ast_GUI

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


global data_empupil; 

data_empupil = [];

% setup parameters
data_empupil.setup.workspace = '';
data_empupil.setup.Pixelsize = 0.12; %um
data_empupil.setup.RefractiveIndex = 1.406;
data_empupil.setup.nMed = 1.352;
data_empupil.setup.Lambda = 0.68;   %um
data_empupil.setup.NA = 1.35;
data_empupil.setup.offset = 100.0;
data_empupil.setup.gain = 2.0;

data_empupil.setup.is_sCMOS = 0;    %0 for EMCCD camera; 1 for sCMOS camera
data_empupil.setup.sCMOS_input = []; 
data_empupil.setup.is_imgsz = 1;    % 0 for image size smaller than 100 x 100 pixels


% pupil parameters
data_empupil.pupil.init_z = zeros(1,21);
data_empupil.pupil.init_z(1) = 1.2;
data_empupil.pupil.ZernikeorderN = 7;  %Zernike order
data_empupil.pupil.Zshift_mode = 1; %XYZ shift mode in each iteration
data_empupil.pupil.Zrange_low = -1;
data_empupil.pupil.Zrange_mid = 0.1;
data_empupil.pupil.Zrange_high = 1;
data_empupil.pupil.Z_pos = [-1:0.1:1];
data_empupil.pupil.bin_lowerBound = 50; %the lower bound of each group images
data_empupil.pupil.iter = 8;
data_empupil.pupil.min_similarity = 0.6; %min similarity in NCC calculation
data_empupil.pupil.blur_sigma = 2;

% reconstruction parameters
data_empupil.recon.isNewdata = 0;   %Reconstruction parameters
data_empupil.recon.isNewpupil = 0;
data_empupil.recon.isSeg = 1;
data_empupil.recon.isRej = 1;
data_empupil.recon.isDC = 1;
data_empupil.recon.isGPU = 0;   %default is CPU version


data_empupil.recon.seg_thresh_low = 25; %segmentation threshold
data_empupil.recon.seg_thresh_high = 40;

data_empupil.recon.rej.min_photon = 1000;   %rejection parameters
data_empupil.recon.rej.llthreshold = 600;
data_empupil.recon.rej.loc_uncer_max = 0.08;
data_empupil.recon.rej.zmask_low = -0.6;
data_empupil.recon.rej.zmask_high = 0.6;

data_empupil.recon.dc.frmpfile = 1000;  %drift correction
data_empupil.recon.dc.z_offset = 2000;
data_empupil.recon.dc.step_ini = 250;

% display parameters
data_empupil.display.Recon_color_highb = 1; %display 
data_empupil.display.imagesz = 256;
data_empupil.display.zm = 5;

%default parameters, load default parameters

% UIWAIT makes INSPR_ast_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = INSPR_ast_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function seg_boxsz_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_boxsz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_boxsz_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_boxsz_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_boxsz_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_boxsz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in data_import_pushbutton.
function data_import_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to data_import_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil
[filename,pathname] = uigetfile(fullfile(data_empupil.setup.workspace,'*.mat'),'Select the data');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.data_show_img_pushbutton,'Enable','on');
   set(handles.data_import_edit,'String', filename);
   set(handles.seg_process_pushbutton,'Enable','on');
   
   % load file
   load([pathname filename]);
   
   if ~exist('ims')
       msgbox('Please check file type: variable should be ims!');
   else
       if isfield(data_empupil,'ims')
           %if import a new data. Make ‘Show subregion’ and ‘Export’ in
           %Segmentation disable. Make 'Show SR' and 'Export' in Display
           %disable
           
           set(handles.seg_export_pushbutton,'Enable','off');
           set(handles.seg_show_img_pushbutton,'Enable','off');
           
           set(handles.display_show_pushbutton,'Enable','off');
           set(handles.display_export_pushbutton,'Enable','off');
       end
       data_empupil.ims = ims(:,:,:);
       data_empupil.display.imagesz = size(ims,1);
       
       if (size(ims,1) > 100)
           data_empupil.setup.is_imgsz = 1;
       else
           data_empupil.setup.is_imgsz = 0;
       end
       msgbox('Finish importing data!');            
   end
end


function data_import_edit_Callback(hObject, eventdata, handles)
% hObject    handle to data_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_import_edit as text
%        str2double(get(hObject,'String')) returns contents of data_import_edit as a double


% --- Executes during object creation, after setting all properties.
function data_import_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in data_show_img_pushbutton.
function data_show_img_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to data_show_img_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

disp('Show projection image of single-molecule dataset...')
% figure; imshow(max(data_empupil.ims,[],3),[]);
figure; imshow(sum(data_empupil.ims,3),[]);
axis tight
title('Sum of diffraction images');



% --- Executes on button press in setup_import_pushbutton.
function setup_import_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setup_import_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

% [filename,pathname] = uigetfile('*.mat','Select the setup configuration');
[filename,pathname] = uigetfile(fullfile(data_empupil.setup.workspace,'*.mat'),'Select the setup configuration');


if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.setup_export_pushbutton,'Enable','on');
   set(handles.setup_import_edit,'String',filename);

   tmp = load([pathname filename]);
   if isfield(tmp,'setup')             
       if isfield(tmp.setup,'Pixelsize')
           data_empupil.setup.Pixelsize = tmp.setup.Pixelsize;
       else
           msgbox('No pixel size. Use the default!');
       end
       
       if isfield(tmp.setup,'RefractiveIndex')
           data_empupil.setup.RefractiveIndex = tmp.setup.RefractiveIndex;
       else
           msgbox('No Refractive Index. Use the default!');
       end
       
       if isfield(tmp.setup,'nMed')
           data_empupil.setup.nMed = tmp.setup.nMed;
       else
           msgbox('No Medium Index. Use the default!');
       end
       
       if isfield(tmp.setup,'Lambda')
           data_empupil.setup.Lambda = tmp.setup.Lambda;
       else
           msgbox('No Lambda. Use the default!');
       end
       
       if isfield(tmp.setup,'NA')
           data_empupil.setup.NA = tmp.setup.NA;
       else
           msgbox('No NA. Use the default!');
       end
       
       if isfield(tmp.setup,'offset')
           data_empupil.setup.offset = tmp.setup.offset;
       else
           msgbox('No camera offset. Use the default!');
       end
       
       if isfield(tmp.setup,'gain')
           data_empupil.setup.gain = tmp.setup.gain;
       else
           msgbox('No camera gain. Use the default!');
       end
       
       msgbox('Finish importing setup configuration!');
   else
       
       msgbox('Please import the correct data!');
   end

end



% --- Executes on button press in setup_set_pushbutton.
function setup_set_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setup_set_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil 
%%
display('Setup parameters...');
% Draw set figure
f_set = figure('Position',[600 100 400 600],'Name','Setup parameters');

%static text
uicontrol(f_set, 'Style','text','String','Pixel size (µm)','FontSize',10,'Position',[50 470 100 50],...
    'HorizontalAlignment','left');
uicontrol(f_set, 'Style','text','String','Refractive index of immersion medium','FontSize',10,'Position',[50 430 115 50],...
    'HorizontalAlignment','left');
uicontrol(f_set, 'Style','text','String','Refractive index of sample medium','FontSize',10,'Position',[50 380 110 50],...
    'HorizontalAlignment','left');
uicontrol(f_set, 'Style','text','String','Lambda (µm)','FontSize',10,'Position',[50 320 100 50],...
    'HorizontalAlignment','left');
uicontrol(f_set, 'Style','text','String','NA','FontSize',10,'Position',[50 270 100 50],...
    'HorizontalAlignment','left');
uicontrol(f_set, 'Style','text','String','Camera offset','FontSize',10,'Position',[50 220 100 50],...
    'HorizontalAlignment','left');
uicontrol(f_set, 'Style','text','String','Camera gain','FontSize',10,'Position',[50 170 100 50],...
    'HorizontalAlignment','left');

%edit 
h_Pixelsize = uicontrol(f_set, 'Style','edit','String',num2str(data_empupil.setup.Pixelsize),...
    'FontSize',10,'Position',[200 500 150 25]);
h_RefractiveIndex = uicontrol(f_set, 'Style','edit','String',num2str(data_empupil.setup.RefractiveIndex),...
    'FontSize',10,'Position',[200 450 150 25]);
h_nMed = uicontrol(f_set, 'Style','edit','String',num2str(data_empupil.setup.nMed),...
    'FontSize',10,'Position',[200 400 150 25]);
h_Lambda  = uicontrol(f_set, 'Style','edit','String',num2str(data_empupil.setup.Lambda),...
    'FontSize',10,'Position',[200 350 150 25]);
h_NA = uicontrol(f_set, 'Style','edit','String',num2str(data_empupil.setup.NA),...
    'FontSize',10,'Position',[200 300 150 25]);
h_offset = uicontrol(f_set, 'Style','edit','String',num2str(data_empupil.setup.offset),...
    'FontSize',10,'Position',[200 250 150 25]);
h_gain = uicontrol(f_set, 'Style','edit','String',num2str(data_empupil.setup.gain),...
    'FontSize',10,'Position',[200 200 150 25]);


% h_zmask_high = uicontrol(f_rej, 'Style','edit','String',num2str(data_empupil.recon.rej.zmask_high,'%0.3f'),...
%    'FontSize',10,'Position',[200 180 70 25]);
%save figure handle
h_set.h_Pixelsize = h_Pixelsize;
h_set.h_RefractiveIndex = h_RefractiveIndex;
h_set.h_nMed = h_nMed;
h_set.h_Lambda = h_Lambda;
h_set.h_NA = h_NA;
h_set.h_offset = h_offset;
h_set.h_gain = h_gain;

% sCMOS part 
h_sCMOS_edit = uicontrol(f_set, 'Style','edit','String','',...
    'FontSize',10,'Position',[200 110 150 25]);

h_sCMOS_import = uicontrol(f_set,'Style','pushbutton','String','Import calibration file','Position',[50 110 130 25],'FontSize',10, 'Callback',...
    {@setup_import_sCMOS_pushbutton_Callback,h_sCMOS_edit});

h_set.h_sCMOS_edit = h_sCMOS_edit;
h_set.h_sCMOS_import = h_sCMOS_import;
h_sCMOS_chk = uicontrol(f_set, 'Style','checkbox','String','sCMOS camera mode','Value',data_empupil.setup.is_sCMOS,...
    'Position',[50 140 200 50],'FontSize',10,...
    'Callback',{@setup_sCMOS_checkBox_Callback,h_set});

h_set.h_sCMOS_chk = h_sCMOS_chk;

if data_empupil.setup.is_sCMOS 
    set(h_set.h_offset,'Enable','off'); %EMCCD off
    set(h_set.h_gain,'Enable','off');
else
    set(h_set.h_sCMOS_edit,'Enable','off'); %sCMOS off
    set(h_set.h_sCMOS_import,'Enable','off');
end

% reset and save button
uicontrol(f_set,'Style','pushbutton','String','Reset','Position',[50 40 100 30],'FontSize',10, 'Callback',...
    {@setup_reset_pushbutton_Callback,h_set});

uicontrol(f_set,'Style','pushbutton','String','Save','Position',[200 40 100 30],'FontSize',10, 'Callback',...
    {@setup_save_pushbutton_Callback,h_set});
%%

set(handles.setup_export_pushbutton,'Enable','on');

%Reset parameters
function setup_sCMOS_checkBox_Callback(src,event,t)
global data_empupil

data_empupil.setup.is_sCMOS = get(src,'Value');

if data_empupil.setup.is_sCMOS
    set(t.h_offset,'Enable','off'); %EMCCD off
    set(t.h_gain,'Enable','off');
    
    set(t.h_sCMOS_edit,'Enable','on'); %sCMOS on
    set(t.h_sCMOS_import,'Enable','on');
else
    set(t.h_offset,'Enable','on');  %EMCCD on
    set(t.h_gain,'Enable','on');
    
    set(t.h_sCMOS_edit,'Enable','off'); %sCMOS off
    set(t.h_sCMOS_import,'Enable','off');
end


function setup_import_sCMOS_pushbutton_Callback(src,event,t)

global data_empupil
[filename,pathname] = uigetfile(fullfile(data_empupil.setup.workspace,'*.mat'),'Select the sCMOS calibration file');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   % load file
   sCMOS_input = load([pathname filename]);   
   set(t,'String', filename);
   
   if isfield(sCMOS_input,'ccdoffset_ch1') && isfield(sCMOS_input,'ccdvar_ch1') && ...
           isfield(sCMOS_input,'gain_ch1')
       data_empupil.setup.sCMOS_input = sCMOS_input;

       msgbox('Finish Import sCMOS calibration file!'); 
   else
       msgbox('sCMOS calibration file import error!');
   end
   
end


%Reset parameters
function setup_reset_pushbutton_Callback(src,event,t)
global data_empupil

default_cfg = load('default_cfg.mat');

%update data_empupil
data_empupil.setup.Pixelsize = default_cfg.setup.Pixelsize;
data_empupil.setup.RefractiveIndex = default_cfg.setup.RefractiveIndex;
data_empupil.setup.nMed = default_cfg.setup.nMed;
data_empupil.setup.Lambda = default_cfg.setup.Lambda;
data_empupil.setup.NA = default_cfg.setup.NA;
data_empupil.setup.offset = default_cfg.setup.offset;
data_empupil.setup.gain = default_cfg.setup.gain;
data_empupil.setup.is_sCMOS = default_cfg.setup.is_sCMOS;
data_empupil.setup.sCMOS_input = default_cfg.setup.sCMOS_input;

%update UI
set(t.h_Pixelsize,'String', num2str(data_empupil.setup.Pixelsize));
set(t.h_RefractiveIndex,'String', num2str(data_empupil.setup.RefractiveIndex));
set(t.h_nMed,'String', num2str(data_empupil.setup.nMed));
set(t.h_Lambda,'String', num2str(data_empupil.setup.Lambda));
set(t.h_NA,'String', num2str(data_empupil.setup.NA));
set(t.h_offset,'String', num2str(data_empupil.setup.offset));
set(t.h_gain,'String', num2str(data_empupil.setup.gain));
set(t.h_offset,'Enable','on'); 
set(t.h_gain,'Enable','on'); 

set(t.h_sCMOS_chk,'Value',data_empupil.setup.is_sCMOS); %sCMOS
set(t.h_sCMOS_import,'Enable','off'); 
set(t.h_sCMOS_edit,'Enable','off'); 
set(t.h_sCMOS_edit,'String', ''); 


%save set parameters
function setup_save_pushbutton_Callback(src,event,t)
global data_empupil 

Pixelsize = str2num( get(t.h_Pixelsize,'String') ); 
RefractiveIndex = str2num( get(t.h_RefractiveIndex,'String') ); 
nMed = str2num( get(t.h_nMed,'String') ); 
Lambda = str2num( get(t.h_Lambda,'String') ); 
NA = str2num( get(t.h_NA,'String') ); 
offset = str2num( get(t.h_offset,'String') ); 
gain = str2num( get(t.h_gain,'String') ); 


data_empupil.setup.Pixelsize = Pixelsize;
data_empupil.setup.RefractiveIndex = RefractiveIndex;
data_empupil.setup.nMed = nMed;
data_empupil.setup.Lambda = Lambda;
data_empupil.setup.NA = NA;
data_empupil.setup.offset = offset;
data_empupil.setup.gain = gain;




function setup_import_edit_Callback(hObject, eventdata, handles)
% hObject    handle to setup_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setup_import_edit as text
%        str2double(get(hObject,'String')) returns contents of setup_import_edit as a double


% --- Executes during object creation, after setting all properties.
function setup_import_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setup_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setup_export_pushbutton.
function setup_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setup_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

[filename, pathname] = uiputfile(fullfile(data_empupil.setup.workspace,'*.mat'), 'save the setup configuration');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   setup = data_empupil.setup;
   save(fullfile(pathname,filename),'setup');
end


function seg_thresh_dist_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_thresh_dist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_thresh_dist_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_thresh_dist_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_thresh_dist_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_thresh_dist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seg_thresh_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_thresh_low_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_thresh_low_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_thresh_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seg_thresh_high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_thresh_high_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_thresh_high_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_thresh_high_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in seg_process_pushbutton.
function seg_process_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to seg_process_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil
addpath('.\Segmentation\');

display('Segmentation process...');
boxsz = str2num( get(handles.seg_boxsz_edit, 'String') );
thresh_dist = str2num( get(handles.seg_thresh_dist_edit, 'String') );
thresh_low = str2num( get(handles.seg_thresh_low_edit, 'String') );
thresh_high = str2num( get(handles.seg_thresh_high_edit, 'String') );

%find subregion
thresh = [thresh_low thresh_high];
[subregion_ch1,seg_display] = crop_subregion_ast(data_empupil.ims,boxsz,thresh,thresh_dist,data_empupil.setup);

data_empupil.subregion_ch1 = subregion_ch1;
data_empupil.seg_display = seg_display;


%show button
num_subregion = size(subregion_ch1,3);
msgbox({'Finish segmentation!' ['There are ' num2str(num_subregion) ' subregion images']});
if num_subregion < 1000
    msgbox('Warning! Less number of sub regions!');
end

display(['Finish segmentation! There are ' num2str(num_subregion) ' subregion images']);

% set(handles.seg_subregion_pushbutton,'Enable','on');
set(handles.seg_export_pushbutton,'Enable','on');
set(handles.seg_show_img_pushbutton,'Enable','on');




% --- Executes on button press in pupil_data_checkbox.
function pupil_data_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_data_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pupil_data_checkbox


if get(handles.pupil_data_checkbox,'Value')==1
   set(handles.pupil_import_pushbutton,'Enable','on');
else
   set(handles.pupil_import_pushbutton,'Enable','off'); 
end


% --- Executes on button press in pupil_import_pushbutton.
function pupil_import_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_import_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data_empupil
[filename,pathname] = uigetfile(fullfile(data_empupil.setup.workspace,'*.mat'),'Select the subregion images..');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.pupil_import_edit,'String',filename);

   % load file
   load([pathname filename]);
    
   if ~exist('subregion_ch1')
       msgbox('Please check sub-region data!');
   else       
       data_empupil.subregion_ch1 = subregion_ch1;
       msgbox('Finish Importing subregion data!');    
   end
   
   
end

function pupil_import_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_import_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_import_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_import_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function pupil_Zrange_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_Zrange_low_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_Zrange_low_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_Zrange_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pupil_Zrange_mid_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_mid_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_Zrange_mid_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_Zrange_mid_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_Zrange_mid_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_mid_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pupil_Zrange_high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_Zrange_high_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_Zrange_high_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_Zrange_high_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pupil_bin_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_bin_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_bin_low_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_bin_low_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_bin_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_bin_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pupil_zshift_mode_popupmenu.
function pupil_zshift_mode_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_zshift_mode_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pupil_zshift_mode_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pupil_zshift_mode_popupmenu
global data_empupil

contents = get(handles.pupil_zshift_mode_popupmenu,'String'); 

switch contents{get(handles.pupil_zshift_mode_popupmenu,'Value')}
    case 'No shift'
        data_empupil.pupil.Zshift_mode = 0;   
    case 'Shift'
        data_empupil.pupil.Zshift_mode = 1;
    otherwise
        data_empupil.pupil.Zshift_mode = 0;
end

display(['Z shift mode is: ' num2str(data_empupil.pupil.Zshift_mode)]);


% --- Executes during object creation, after setting all properties.
function pupil_zshift_mode_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_zshift_mode_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_data_checkbox.
function recon_data_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_data_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_data_checkbox
global data_empupil

if get(handles.recon_data_checkbox,'Value')==1
    data_empupil.recon.isNewdata = 1;
    set(handles.recon_import_data_pushbutton,'Enable','on');
else
    data_empupil.recon.isNewdata = 0; 
    set(handles.recon_import_data_pushbutton,'Enable','off');
end


% --- Executes on button press in recon_import_data_pushbutton.
function recon_import_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data_empupil
[filename,pathname] = uigetfile(fullfile(data_empupil.setup.workspace,'*.mat'),'Select the data (you can choose multiple data)...','MultiSelect','on');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.recon_import_data_edit,'String', pathname);
   
   % load one file or muti-files
   if iscell(filename)
       dirN = numel(filename);
       data_empupil.recon.datafile_name = filename;
   else
       dirN = 1;
       data_empupil.recon.datafile_name{1} = filename;
   end
   
   %save file path
   data_empupil.recon.dirN = dirN;
   data_empupil.recon.datapath = pathname;
   
   msgbox({ 'Finish Importing the data path!' ['You choose ' num2str(dirN) ' data'] });

   
   
end

function recon_import_data_edit_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_import_data_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_import_data_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_import_data_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_import_data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function recon_import_tform_edit_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_tform_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_import_tform_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_import_tform_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_import_tform_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_import_tform_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_pupil_checkbox.
function recon_pupil_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_pupil_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_pupil_checkbox
global data_empupil

if get(handles.recon_pupil_checkbox,'Value')==1
    data_empupil.recon.isNewpupil = 1;
    set(handles.recon_import_pupil_pushbutton,'Enable','on');
else
    data_empupil.recon.isNewpupil = 0;

    set(handles.recon_import_pupil_pushbutton,'Enable','off');
end



% --- Executes on button press in recon_import_pupil_pushbutton.
function recon_import_pupil_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_pupil_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil
[filename,pathname] = uigetfile(fullfile(data_empupil.setup.workspace,'*.mat'),'Select the pupil model (you can choose multiple pupil)...','MultiSelect','on');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.recon_import_pupil_edit,'String', pathname);
   
   % load one file or muti-files
   probj_all = [];
   loopn = [];
   if iscell(filename)
       loopn = numel(filename);
       for ii = 1 : loopn
           input_file = [pathname filename{ii}];
           tmp = load(input_file);
           if isfield(tmp,'probj')
               probj_all{ii} = tmp.probj;
           else
               msgbox('Please import the correct data!');
               return;
           end
       end
       msgbox({'Finish importing the pupil model!' ['You choose ' num2str(loopn) ' pupil model']});
   else
       loopn = 1;
       input_file = [pathname filename];
       tmp = load(input_file);
       if isfield(tmp,'probj')                 
           probj_all{1} = tmp.probj;
       else
           msgbox('Please import the correct data!');
           return;
       end
       msgbox({'Finish importing the pupil model!' ['You choose ' num2str(loopn) ' pupil model']});
   end
   
   %save all pupil models to data_empupil
   data_empupil.recon.probj_all = probj_all;
   data_empupil.recon.loopn = loopn;
   
end






function recon_import_pupil_edit_Callback(hObject, ~, handles)
% hObject    handle to recon_import_pupil_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_import_pupil_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_import_pupil_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_import_pupil_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_import_pupil_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_seg_checkbox.
function recon_seg_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_seg_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_seg_checkbox
global data_empupil

if get(handles.recon_seg_checkbox,'Value')==1
    set(handles.recon_seg_thresh_low_edit,'Enable','on');
    set(handles.recon_seg_thresh_high_edit,'Enable','on');
    
    data_empupil.recon.isSeg = 1;
else
    set(handles.recon_seg_thresh_low_edit,'String',num2str(data_empupil.recon.seg_thresh_low));
    set(handles.recon_seg_thresh_high_edit,'String',num2str(data_empupil.recon.seg_thresh_high));
    set(handles.recon_seg_thresh_low_edit,'Enable','off');
    set(handles.recon_seg_thresh_high_edit,'Enable','off');
    
    data_empupil.recon.isSeg = 0;
end




function recon_seg_thresh_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to recon_seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_seg_thresh_low_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_seg_thresh_low_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_seg_thresh_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function recon_seg_thresh_high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to recon_seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_seg_thresh_high_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_seg_thresh_high_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_seg_thresh_high_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_rej_checkbox.
function recon_rej_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_rej_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_rej_checkbox
global data_empupil

if get(handles.recon_rej_checkbox,'Value')==1
    data_empupil.recon.isRej = 1;
    
    set(handles.recon_rej_change_pushbutton,'Enable','on');
else
    data_empupil.recon.isRej = 0;
    
    data_empupil.recon.rej.min_photon = 1000;
    data_empupil.recon.rej.llthreshold = 1000;
    data_empupil.recon.rej.loc_uncer_max = 0.08;
    data_empupil.recon.rej.zmask_low = -0.8;
    data_empupil.recon.rej.zmask_high = 0.8;
    set(handles.recon_rej_change_pushbutton,'Enable','off');
end


% --- Executes on button press in recon_rej_change_pushbutton.
function recon_rej_change_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_rej_change_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil 

% Draw rejection figure
f_rej = figure('Position',[600 200 400 400],'Name','Rejection parameters');

%static text
uicontrol(f_rej, 'Style','text','String','Min photon','FontSize',10,'Position',[50 300 100 50],...
    'HorizontalAlignment','left');
uicontrol(f_rej, 'Style','text','String','LLR','FontSize',10,'Position',[50 250 100 50],...
    'HorizontalAlignment','left');
uicontrol(f_rej, 'Style','text','String','Max Z uncertainty (µm)','FontSize',10,'Position',[50 200 100 50],...
    'HorizontalAlignment','left');
uicontrol(f_rej, 'Style','text','String','Z mask (µm)','FontSize',10,'Position',[50 150 100 50],...
    'HorizontalAlignment','left');

%edit 
h_photon  = uicontrol(f_rej, 'Style','edit','String',num2str(data_empupil.recon.rej.min_photon),...
    'FontSize',10,'Position',[180 330 150 25]);
h_ll = uicontrol(f_rej, 'Style','edit','String',num2str(data_empupil.recon.rej.llthreshold),...
    'FontSize',10,'Position',[180 280 150 25]);
h_uncer = uicontrol(f_rej, 'Style','edit','String',num2str(data_empupil.recon.rej.loc_uncer_max,'%0.3f'),...
    'FontSize',10,'Position',[180 230 150 25]);
h_zmask_low = uicontrol(f_rej, 'Style','edit','String',num2str(data_empupil.recon.rej.zmask_low,'%0.3f'),...
    'FontSize',10,'Position',[180 180 70 25]);
h_zmask_high = uicontrol(f_rej, 'Style','edit','String',num2str(data_empupil.recon.rej.zmask_high,'%0.3f'),...
    'FontSize',10,'Position',[260 180 70 25]);

h_rej.h_photon = h_photon;
h_rej.h_ll = h_ll;
h_rej.h_uncer = h_uncer;
h_rej.h_zmask_low = h_zmask_low;
h_rej.h_zmask_high = h_zmask_high;

%save button
uicontrol(f_rej,'Style','pushbutton','String','Save','Position',[180 80 100 30],'FontSize',10,'Callback',...
    {@recon_rej_save_pushbutton_Callback,h_rej});

%Reset button
uicontrol(f_rej,'Style','pushbutton','String','Reset','Position',[50 80 100 30],'FontSize',10,'Callback',...
    {@recon_rej_reset_pushbutton_Callback,h_rej});


%Reset parameters
function recon_rej_reset_pushbutton_Callback(src,event,t)
global data_empupil

default_cfg = load('default_cfg.mat');

%update data_empupil
data_empupil.recon.rej.min_photon = default_cfg.recon.rej.min_photon;
data_empupil.recon.rej.llthreshold = default_cfg.recon.rej.llthreshold;
data_empupil.recon.rej.loc_uncer_max = default_cfg.recon.rej.loc_uncer_max;
data_empupil.recon.rej.zmask_low = default_cfg.recon.rej.zmask_low;
data_empupil.recon.rej.zmask_high = default_cfg.recon.rej.zmask_high;


%update UI
set(t.h_photon,'String', num2str(data_empupil.recon.rej.min_photon));
set(t.h_ll,'String', num2str(data_empupil.recon.rej.llthreshold));
set(t.h_uncer,'String', num2str(data_empupil.recon.rej.loc_uncer_max));
set(t.h_zmask_low,'String', num2str(data_empupil.recon.rej.zmask_low));
set(t.h_zmask_high,'String', num2str(data_empupil.recon.rej.zmask_high));


%save rejection parameters
function recon_rej_save_pushbutton_Callback(src,event,t)
global data_empupil 

min_photon = str2num( get(t.h_photon,'String') ); 
llthreshold = str2num( get(t.h_ll,'String') ); 
loc_uncer_max = str2num( get(t.h_uncer,'String') ); 
zmask_low = str2num( get(t.h_zmask_low,'String') ); 
zmask_high = str2num( get(t.h_zmask_high,'String') ); 

data_empupil.recon.rej.min_photon = min_photon;
data_empupil.recon.rej.llthreshold = llthreshold;
data_empupil.recon.rej.loc_uncer_max = loc_uncer_max;
data_empupil.recon.rej.zmask_low = zmask_low;
data_empupil.recon.rej.zmask_high = zmask_high;





% --- Executes on button press in recon_dc_checkbox.
function recon_dc_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_dc_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_dc_checkbox

global data_empupil

if get(handles.recon_dc_checkbox,'Value')==1
    data_empupil.recon.isDC = 1;
    
    set(handles.recon_dc_change_pushbutton,'Enable','on');
else
    data_empupil.recon.isDC = 0;
    
    data_empupil.recon.dc.frmpfile = 1000;
    data_empupil.recon.dc.z_offset = 2000;
    data_empupil.recon.dc.step_ini = 400;
    set(handles.recon_dc_change_pushbutton,'Enable','off');
end




% --- Executes on button press in recon_dc_change_pushbutton.
function recon_dc_change_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_dc_change_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil 

% Draw DC figure
f_dc = figure('Position',[600 200 350 300],'Name','Drift correction parameters');

% static text
uicontrol(f_dc, 'Style','text','String','Frame bin','FontSize',10,'Position',[50 200 100 50],...
    'HorizontalAlignment','left');
uicontrol(f_dc, 'Style','text','String','Initial Z offset (nm)','FontSize',10,'Position',[50 150 150 50],...
    'HorizontalAlignment','left');
uicontrol(f_dc, 'Style','text','String','Step interval (nm)','FontSize',10,'Position',[50 100 150 50],...
    'HorizontalAlignment','left');

% edit 
h_frmp  = uicontrol(f_dc, 'Style','edit','String',num2str(data_empupil.recon.dc.frmpfile),...
    'FontSize',10,'Position',[180 230 110 25]);
h_z_offset = uicontrol(f_dc, 'Style','edit','String',num2str(data_empupil.recon.dc.z_offset),...
    'FontSize',10,'Position',[180 180 110 25]);
h_step = uicontrol(f_dc, 'Style','edit','String',num2str(data_empupil.recon.dc.step_ini),...
    'FontSize',10,'Position',[180 130 110 25]);

% data_empupil.recon.dc.frmpfile
% data_empupil.recon.dc.z_offset
% data_empupil.recon.dc.step_ini

h_dc.h_frmp = h_frmp;
h_dc.h_z_offset = h_z_offset;
h_dc.h_step = h_step;

% save button
uicontrol(f_dc,'Style','pushbutton','String','Save','Position',[180 50 100 30],'FontSize',10,'Callback',...
    {@recon_dc_save_pushbutton_Callback,h_dc});
% Reset button
uicontrol(f_dc,'Style','pushbutton','String','Reset','Position',[50 50 100 30],'FontSize',10,'Callback',...
    {@recon_dc_reset_pushbutton_Callback,h_dc});

%Reset parameters
function recon_dc_reset_pushbutton_Callback(src,event,t)
global data_empupil

default_cfg = load('default_cfg.mat');

%update data_empupil
data_empupil.recon.dc.frmpfile = default_cfg.recon.dc.frmpfile;
data_empupil.recon.dc.z_offset = default_cfg.recon.dc.z_offset;
data_empupil.recon.dc.step_ini = default_cfg.recon.dc.step_ini;


%update UI
set(t.h_frmp,'String', num2str(data_empupil.recon.dc.frmpfile));
set(t.h_z_offset,'String', num2str(data_empupil.recon.dc.z_offset));
set(t.h_step,'String', num2str(data_empupil.recon.dc.step_ini));


%save rejection parameters
function recon_dc_save_pushbutton_Callback(src,event,t)
global data_empupil 

frmpfile = str2num( get(t.h_frmp,'String') ); 
z_offset = str2num( get(t.h_z_offset,'String') ); 
step_ini = str2num( get(t.h_step,'String') ); 

data_empupil.recon.dc.frmpfile = frmpfile;
data_empupil.recon.dc.z_offset = z_offset;
data_empupil.recon.dc.step_ini = step_ini;


% --- Executes on button press in recon_process_pushbutton.
function recon_process_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_process_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('.\3D_localization\');

global data_empupil
global recon_stop

recon_stop = 0; %initil stop label is none

% check the data
if data_empupil.recon.isNewdata == 0
    if ~isfield(data_empupil,'ims')
        msgbox('Please import the data!');
        return;
    end
    data_empupil.recon.ims = data_empupil.ims;
    data_empupil.recon.dirN = 1;
    data_empupil.recon.datapath = pwd;
else
    if ~isfield(data_empupil.recon,'datapath')
        msgbox('Please import the data!');
        return;
    end
end


if data_empupil.recon.isNewpupil == 0
    if ~isfield(data_empupil,'probj')
        msgbox('Please import pupil model!');
        return;
    end
    data_empupil.recon.probj_all{1} = data_empupil.probj;
    data_empupil.recon.loopn = 1;
else
    if ~isfield(data_empupil.recon,'probj_all')
        msgbox('Please import the pupil model!');
        return;
    end
end

if data_empupil.recon.isSeg == 1
    data_empupil.recon.seg_thresh_low = str2num( get(handles.recon_seg_thresh_low_edit, 'String') );
    data_empupil.recon.seg_thresh_high = str2num( get(handles.recon_seg_thresh_high_edit, 'String') );    
end

if data_empupil.setup.is_sCMOS
    disp('Consider independent readout noise in sCMOS camera');
end

if data_empupil.recon.isGPU
    %CUDA check
    try, listGPUs;
    catch
        msgbox('Please install CUDA environment! https://developer.nvidia.com/cuda-75-downloads-archive');
        return;
    end  
end

%SMLM pupil fitting
srobj = analysis3D_fromPupil_ast(data_empupil.recon,data_empupil.setup);

if recon_stop == 1
    msgbox('Stop by user control!');
    return;
end

data_empupil.srobj = srobj; 

% set export
msgbox('Finish SMLM reconstruction!');
set(handles.recon_export_pushbutton,'Enable','on');
set(handles.recon_export_csv_pushbutton,'Enable','on');

% make dispaly buttons enable
set(handles.display_show_pushbutton,'Enable','on');
set(handles.display_export_pushbutton,'Enable','off');




% --- Executes on button press in display_set_pushbutton.
function display_set_pushbutton_Callback(~, eventdata, handles)
% hObject    handle to display_set_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

% Draw DC figure
f_display = figure('Position',[600 200 350 260],'Name','Set display parameters');

% static text
uicontrol(f_display, 'Style','text','String','Image saturation','FontSize',10,'Position',[50 150 100 50],...
    'HorizontalAlignment','left');
% uicontrol(f_display, 'Style','text','String','Image width (pixels)','FontSize',10,'Position',[50 150 150 50],...
%     'HorizontalAlignment','left');
uicontrol(f_display, 'Style','text','String','Zoom in','FontSize',10,'Position',[50 100 100 50],...
    'HorizontalAlignment','left');

% edit 
h_contrast  = uicontrol(f_display, 'Style','edit','String',num2str(data_empupil.display.Recon_color_highb),...
    'FontSize',10,'Position',[180 180 110 25]);
% h_imgsz = uicontrol(f_display, 'Style','edit','String',num2str(data_empupil.display.imagesz),...
%     'FontSize',10,'Position',[180 180 110 25]);
h_zm = uicontrol(f_display, 'Style','edit','String',num2str(data_empupil.display.zm),...
    'FontSize',10,'Position',[180 130 110 25]);

h_display.h_contrast = h_contrast;
% h_display.h_imgsz = h_imgsz;
h_display.h_zm = h_zm;

% save button
uicontrol(f_display,'Style','pushbutton','String','save','Position',[180 50 100 30],'FontSize',10,'Callback',...
    {@display_set_save_pushbutton_Callback,h_display});

%save display parameters
function display_set_save_pushbutton_Callback(src,event,t)
global data_empupil 

Recon_color_highb = str2num( get(t.h_contrast,'String') ); 
% imagesz = str2num( get(t.h_imgsz,'String') ); 
zm = str2num( get(t.h_zm,'String') ); 

data_empupil.display.Recon_color_highb = Recon_color_highb;
% data_empupil.display.imagesz = imagesz;
data_empupil.display.zm = zm;


% --- Executes on button press in display_show_pushbutton.
function display_show_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to display_show_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil 

if ~isfield(data_empupil,'srobj')
    msgbox('No reconstruction results!');
    return;
end

display('Show 3D reconstrcted Z-color image...');

obj = data_empupil.srobj;
sz = data_empupil.display.imagesz;
obj.zm = data_empupil.display.zm;
obj.Recon_color_highb = data_empupil.display.Recon_color_highb;
segnum = 64;
flagstr = [];
if isempty(obj.loc_x_f)||isempty(obj.loc_y_f)||isempty(obj.loc_z_f)
    warning('loc_x(y or z)_f property is empty, reconstructing from raw localization data');
    reconx=obj.loc_x;
    recony=obj.loc_y;
    reconz=obj.loc_z;
    flagstr='raw';
else
    reconx=obj.loc_x_f(:)./obj.Cam_pixelsz;
    recony=obj.loc_y_f(:)./obj.Cam_pixelsz;
    reconz=obj.loc_z_f(:);
    flagstr='dc';
end

if isempty(reconx)||isempty(recony)||isempty(reconz)
    error('Empty matrix detected. Reconstruction will not proceed.');
end

[rch,gch,bch]=srhist_color(sz,obj.zm,reconx,recony,reconz,segnum);
% save colored reconstruction
rchsm = imgaussfilt(rch,1);
gchsm = imgaussfilt(gch,1);
bchsm = imgaussfilt(bch,1);
rchsmst=imstretch_linear(rchsm,0,obj.Recon_color_highb,0,255);
gchsmst=imstretch_linear(gchsm,0,obj.Recon_color_highb,0,255);
bchsmst=imstretch_linear(bchsm,0,obj.Recon_color_highb,0,255);
colorim = cat(3,rchsmst,gchsmst,bchsmst);
colorim = uint8(colorim);

figure; imshow(colorim,[]);
axis tight

data_empupil.display.colorim = colorim;
% save backup
imwrite(colorim,fullfile(data_empupil.setup.workspace,['SR_Zcol_' flagstr '_backup.tif']));


%enable export button
set(handles.display_export_pushbutton,'Enable','on');



% --- Executes on button press in recon_export_csv_pushbutton.
function recon_export_csv_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_export_csv_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

if ~isfield(data_empupil,'srobj')
    msgbox('No reconstruction results!');
    return;
end

obj = data_empupil.srobj;
[filename, pathname] = uiputfile(fullfile(data_empupil.setup.workspace,'particles.csv'), 'Export to CSV format...');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])
   
%    flagstr=[];
   if isempty(obj.loc_x_f)||isempty(obj.loc_y_f)||isempty(obj.loc_z_f)
       warning('loc_x(y or z)_f property is empty, exporting from raw localization data');
       reconx=obj.loc_x;
       recony=obj.loc_y;
       reconz=obj.loc_z;
%        flagstr='raw';
   else
       reconx=obj.loc_x_f./obj.Cam_pixelsz;
       recony=obj.loc_y_f./obj.Cam_pixelsz;
       reconz=obj.loc_z_f;
%        flagstr='dc';
   end
   
   [flag]=export2csv(pathname,filename,1,{reconx.*obj.Cam_pixelsz},{recony.*obj.Cam_pixelsz},{reconz},{obj.loc_t});
      
   if flag==1
       msgbox('Export Successfully!');
   else
       msgbox('Export with error!');
   end

end




function pupil_iter_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_iter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_iter_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_iter_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_iter_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_iter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pupil_min_similarity_edit_Callback(hObject, eventdata, ~)
% hObject    handle to pupil_min_similarity_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_min_similarity_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_min_similarity_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_min_similarity_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_min_similarity_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pupil_zernike_order_popupmenu.
function pupil_zernike_order_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_zernike_order_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pupil_zernike_order_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pupil_zernike_order_popupmenu
global data_empupil

contents = get(handles.pupil_zernike_order_popupmenu,'String'); 

switch contents{get(handles.pupil_zernike_order_popupmenu,'Value')}
    case '25'
        data_empupil.pupil.ZernikeorderN = 4;   
    case '36'
        data_empupil.pupil.ZernikeorderN = 5;
    case '49'
        data_empupil.pupil.ZernikeorderN = 6;
    case '64'
        data_empupil.pupil.ZernikeorderN = 7;
    otherwise
        data_empupil.pupil.ZernikeorderN = 7;
end

display(['shift mode is: ' num2str(data_empupil.pupil.ZernikeorderN)]);


% --- Executes during object creation, after setting all properties.
function pupil_zernike_order_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_zernike_order_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pupil_change_pushbutton.
function pupil_change_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_change_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

f_init = figure('Position',[600 400 1000 150],'Name','Initial Zernike value');

cnames = {'Ast','DAst','Coma x','Coma y','1st Sph','Trefoil x','Trefoil y','2nd Ast','2nd DAst',...
    '2nd Coma x','2nd Coma y','2nd Sph','Tetrafoil x','Tetrafoil y','2nd Trefoil x','2nd Trefoil y',...
    '3rd Ast','3rd DAst','3rd Coma x','3rd Coma y','3rd Sph'};
rnames = {'Value'};
dat_zernike = data_empupil.pupil.init_z; 

t_zernike = uitable(f_init,'Data',dat_zernike,'ColumnName',cnames,'RowName',rnames,'Position',[20 20 800 100]);
t_zernike.ColumnEditable = true;

% data_empupil.pupil.init_z

uicontrol(f_init,'Style','pushbutton','String','save','Position',[860 50 100 50],'Callback',{@pupil_zernike_save_pushbutton_Callback,t_zernike});

%save initial zernike value
function pupil_zernike_save_pushbutton_Callback(src,event,t)
global data_empupil 

data_empupil.pupil.init_z = get(t,'Data');
display(['initial zernike is: ' num2str(data_empupil.pupil.init_z)]);



function pupil_blur_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_blur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_blur_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_blur_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_blur_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_blur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pupil_process_pushbutton.
function pupil_process_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_process_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('.\INSPR_model_generation\');

global data_empupil
global pupil_stop

%check data

if ~isfield(data_empupil,'subregion_ch1')
    msgbox('Please improt subregion images!');
    return
end

pupil_stop = 0; %initial stop == 0


%import parameters from GUI
display('Import parameters...');

Zrange_low = str2num(get(handles.pupil_Zrange_low_edit,'String'));
Zrange_mid = str2num(get(handles.pupil_Zrange_mid_edit,'String'));
Zrange_high = str2num(get(handles.pupil_Zrange_high_edit,'String'));
bin_lowerBound = str2num(get(handles.pupil_bin_low_edit,'String'));
iter = str2num(get(handles.pupil_iter_edit,'String'));
min_similarity = str2num(get(handles.pupil_min_similarity_edit,'String'));
% blur_sigma = str2num(get(handles.pupil_blur_edit,'String'));

Z_pos = [Zrange_low:Zrange_mid:Zrange_high];

%import to data_empupil
data_empupil.pupil.Zrange_low = Zrange_low;
data_empupil.pupil.Zrange_mid = Zrange_mid;
data_empupil.pupil.Zrange_high = Zrange_high;
data_empupil.pupil.Z_pos = Z_pos;
data_empupil.pupil.bin_lowerBound = bin_lowerBound; %the lower bound of each group images
data_empupil.pupil.iter = iter;
data_empupil.pupil.min_similarity = min_similarity; %min similarity in NCC calculation

%normalization for subregion images
display('Image normalization...');

subregion_ch1_norm = subregion_normalization_ast(data_empupil.subregion_ch1);

%EMpupil process... unfinished
display('EMpupil process...');


probj = INSPR_model_generation_ast(subregion_ch1_norm,data_empupil.setup,data_empupil.pupil); %this function is initial, need change

if pupil_stop == 1
    msgbox('Stop by user control!'); 
    return;
end

%save probj (pupil)
data_empupil.probj = probj;


%enable button
msgbox('Finish pupil estimation!');

set(handles.pupil_export_pushbutton,'Enable','on');
set(handles.pupil_showPSFs_pushbutton,'Enable','on');
set(handles.pupil_show_pupil_pushbutton,'Enable','on');
set(handles.pupil_show_zernike_pushbutton,'Enable','on');



% --- Executes on button press in pupil_export_pushbutton.
function pupil_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data_empupil

[filename, pathname] = uiputfile(fullfile(data_empupil.setup.workspace,'*.mat'), 'save pupil...');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   probj = data_empupil.probj;
   save(fullfile(pathname,filename),'probj');
end


% --- Executes on button press in pupil_showPSFs_pushbutton.
function pupil_showPSFs_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_showPSFs_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

disp('Show the Reassembled and EMpupil PSFs ...')
% data_empupil.probj.genPRfigs('PSF');
genPupilfigs(data_empupil.probj, 'PSF',data_empupil.setup.workspace);

% --- Executes on button press in pupil_show_pupil_pushbutton.
function pupil_show_pupil_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_show_pupil_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

disp('Show PR pupil and Zernike pupil ...')
% data_empupil.probj.genPRfigs('pupil');
genPupilfigs(data_empupil.probj, 'pupil',data_empupil.setup.workspace);

% --- Executes on button press in pupil_show_zernike_pushbutton.
function pupil_show_zernike_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_show_zernike_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data_empupil

disp('Show Zernike value in magnitude and phase ...')
% data_empupil.probj.genPRfigs('zernike');
genPupilfigs(data_empupil.probj, 'zernike',data_empupil.setup.workspace);


% --- Executes on button press in seg_export_pushbutton.
function seg_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to seg_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data_empupil

[filename, pathname] = uiputfile(fullfile(data_empupil.setup.workspace,'*.mat'), 'save subregion images');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   subregion_ch1 = data_empupil.subregion_ch1;

   save(fullfile(pathname,filename),'subregion_ch1');
end



% --- Executes on button press in pupil_zernike_checkbox.
function pupil_zernike_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_zernike_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pupil_zernike_checkbox
global data_empupil

if get(handles.pupil_zernike_checkbox,'Value')==1
    set(handles.pupil_change_pushbutton,'Enable','on');
else
    data_empupil.pupil.init_z = zeros(1,21);    %initial pupil with prior knowledge
    data_empupil.pupil.init_z(1) = 1.2;
    set(handles.pupil_change_pushbutton,'Enable','off');
end




% --- Executes on button press in recon_export_pushbutton.
function recon_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data_empupil

[filename, pathname] = uiputfile(fullfile(data_empupil.setup.workspace,'*.mat'), 'save SMLM reconstruction result...');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   srobj = data_empupil.srobj;
   save(fullfile(pathname,filename),'srobj');
   msgbox('Finish export!');
end


% --- Executes on button press in pupil_stop_pushbutton.
function pupil_stop_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_stop_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pupil_stop

pupil_stop = 1;


% --- Executes on button press in recon_stop_pushbutton.
function recon_stop_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_stop_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global recon_stop

recon_stop = 1;


% --- Executes on button press in setup_workspace_pushbutton.
function setup_workspace_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setup_workspace_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data_empupil

folder_name = uigetdir(data_empupil.setup.workspace,'Select workspace folder');

if isequal(folder_name,0)
   disp('User selected Cancel')
else
   disp(['User selected ', folder_name]);
   set(handles.setup_workspace_edit,'String', folder_name);

   data_empupil.setup.workspace = folder_name;
end




function setup_workspace_edit_Callback(hObject, eventdata, handles)
% hObject    handle to setup_workspace_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setup_workspace_edit as text
%        str2double(get(hObject,'String')) returns contents of setup_workspace_edit as a double


% --- Executes during object creation, after setting all properties.
function setup_workspace_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setup_workspace_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in display_export_pushbutton.
function display_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to display_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

if isfield(data_empupil.display,'colorim')
    [filename, pathname] = uiputfile(fullfile(data_empupil.setup.workspace,'*.tif'), ...
        'save super resolution image...');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User selected Cancel')
    else
        disp(['User selected ',fullfile(pathname,filename)])
                
        imwrite(data_empupil.display.colorim,fullfile(pathname,filename),'TIFF');
        msgbox('Finish export!');
    end

else
    msgbox('No Z projection super-resolution images!');
end

% data_empupil.display.colorim = colorim;



function seg_frame_id_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_frame_id_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_frame_id_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_frame_id_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_frame_id_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_frame_id_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in seg_show_img_pushbutton.
function seg_show_img_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to seg_show_img_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil

boxsz = str2num( get(handles.seg_boxsz_edit, 'String') );
num_display = str2num( get(handles.seg_frame_id_edit, 'String') );

if num_display < 1 || num_display > size(data_empupil.seg_display.ims_ch1,3)
    msgbox('Please input the correct number!');
end

display('Show segmentation results');
raw = data_empupil.seg_display.ims_ch1(:,:,num_display)/max(max(data_empupil.seg_display.ims_ch1(:,:,num_display)));
index_rec = find(data_empupil.seg_display.allcds_mask(:,3) == num_display-1);

rec_vector = cat(2,data_empupil.seg_display.t1(index_rec),data_empupil.seg_display.l1(index_rec),repmat(boxsz,length(index_rec),2));
img_select = insertShape(raw,'Rectangle',rec_vector,'LineWidth',1, 'Color', 'green');

figure; imshow(img_select);
axis tight
title('Segmentation results');



% --- Executes on button press in recon_gpu_checkbox.
function recon_gpu_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_gpu_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_gpu_checkbox
global data_empupil

if get(handles.recon_gpu_checkbox,'Value')==1
    data_empupil.recon.isGPU = 1;
else
    data_empupil.recon.isGPU = 0;
end


% --- Executes on button press in display_import_pushbutton.
function display_import_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to display_import_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_empupil
[filename,pathname] = uigetfile(fullfile(data_empupil.setup.workspace,'*.mat'),'Select reconstruction result');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   % load file
   tmp = load([pathname filename]);
   
   if isfield(tmp,'srobj')        
       data_empupil.srobj = tmp.srobj;
       set(handles.display_show_pushbutton,'Enable','on');
       set(handles.display_export_pushbutton,'Enable','off');
       set(handles.display_import_edit,'String', filename);
       
       msgbox('Finish importing reconstruction result!');
       
   else
       msgbox('Please import the correct data!');
   end
end


function display_import_edit_Callback(hObject, eventdata, handles)
% hObject    handle to display_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of display_import_edit as text
%        str2double(get(hObject,'String')) returns contents of display_import_edit as a double


% --- Executes during object creation, after setting all properties.
function display_import_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
