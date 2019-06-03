function varargout = Set_Mission(varargin)
% SET_MISSION MATLAB code for Set_Mission.fig
%      SET_MISSION, by itself, creates a new SET_MISSION or raises the existing
%      singleton*.
%
%      H = SET_MISSION returns the handle to a new SET_MISSION or the handle to
%      the existing singleton*.
%
%      SET_MISSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SET_MISSION.M with the given input arguments.
%
%      SET_MISSION('Property','Value',...) creates a new SET_MISSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Set_Mission_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Set_Mission_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Set_Mission

% Last Modified by GUIDE v2.5 22-Sep-2017 12:47:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Set_Mission_OpeningFcn, ...
                   'gui_OutputFcn',  @Set_Mission_OutputFcn, ...
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


% --- Executes just before Set_Mission is made visible.
function Set_Mission_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Set_Mission (see VARARGIN)

% Choose default command line output for Set_Mission
handles.output = hObject;

global This;    % Project Object currently being used

global mint;    % Handle for Minimum Time Listbox 
global smint;   % String array for Minimum Time Listbox

global maxt;    % Handle for Maximum Time Listbox
global smaxt;   % Sring array for Maximum Time Listbox

global nomt;    % Handle for Nominal Time (Initial Guess) Listbox
global snomt;   % Sring array for Nominal Time Listbox

global perh;    % Handle for Minimum Perihelion Distance Listbox
global sperhel; % String array for Minimum Perihelion Listbox

global per;     % Handle for Encounter Constraint Details Minimum Periapsis altitude Listbox
global per2;    % Handle for Encounter Constraint Details Maximum DeltaV Listbox
global sper;    % 2-Dim String array for Encounter Constraint Details Listboxes

global plan;    % Handle for Planet name Listbox
global splan;   % String array for Planet name Listbox

global dur;     % Handle for Maximum Duration Edit box
global durstring;   % String for Maximum Duration Edit box

global run;     % Handle for Maximum Optimization Time Limit Edit box


%
% Firstly Construct Edit Box for Max Optimization Time
%
run = uicontrol('Style', 'edit','Callback',@Run_Time_Set);
run.FontSize=14.0;
run.String = sprintf("%5.2f", This.Run_Time/60 );
run.Position = [550 90 100 40];

%
% Prepare String for Max Mission Duration Edit Box
%
if (This.Max_Duration >= This.MAX_DURATION)
    durstring = "MAX";
else
    durstring = sprintf("%5.2f", This.Max_Duration/365/24/60/60);
end

%
% Now Construct Edit Box for Max Mission Duration
%
dur = uicontrol('Style', 'edit','Callback',@Max_Dur_Set);
dur.FontSize=14.0;
dur.String = durstring;
dur.Position = [ 200 90 100 40];

%
% Now Construct flyby/Rendezvous Checkbox
%
flyren = uicontrol('Style', 'checkbox', 'Value',This.Current_Mission.FlybyRendez,'Callback', @Destination_Arrival_Set);
    flyren.Position=[900 90 100 40];
    flyren.FontSize=25.0;

%
% Now Construct Prograde or otherwise Checkbox
%    
prograde = uicontrol('Style', 'checkbox', 'Value',This.Current_Mission.wayflag,'Callback', @Prograde_Set);
    prograde.Position=[1250 90 100 40];
    prograde.FontSize=35.0;

%
% Prepare String for Perihelion Listbox
%
sperhel=strings(This.Current_Mission.Trajectory.Nbody);

for i=1:This.Current_Mission.Trajectory.Nbody
    sperhel(i)=sprintf("%12.6f",This.Perihelia(i)/This.AU);
    if (bitand(This.Perihelia_flag,2^(i-1)))
        sperhel(i)='P'+sperhel(i);
    end
end

%
% Now construct Maximum Perihelion Listbox
%
perh=uicontrol('Style', 'listbox','Callback',@Perihelia_Set);
perh.FontSize=12;
perh.String=sperhel;
perh.Position=[400 230 120 290]; 

%
% Prepare String for Planet name Listbox
%
splan=strings(This.Current_Mission.Trajectory.Nbody);

for i=1:This.Current_Mission.Trajectory.Nbody
    splan(i)=This.Current_Mission.Trajectory.Body_Set(i).name;
end

%
% Now construct Planet name Listbox
%
plan=uicontrol('Style', 'listbox','Callback',@Planet_Name_Set);
plan.FontSize=12;
plan.String=splan;
plan.Position=[90 230 120 290];

%
% Prepare string for Solar System Object SPICE ID Listbox
%
sID=strings(This.Current_Mission.Trajectory.Nbody);

for i=1:This.Current_Mission.Trajectory.Nbody
    sID(i)=This.Current_Mission.Trajectory.Body_Set(i).ID;
end

%
% Now Construct handle for SPICE ID Listbox NOTE ID IS NOT MODIFIABLE BY
% USER
%
ID=uicontrol('Style', 'listbox');
ID.FontSize=12;
ID.String=sID;
ID.Position=[230 230 120 290];

%
% Prepare String array for Minimum Allowable Mission Times for Optimization
%
smint=strings(This.Current_Mission.Trajectory.Nbody);

for i=1:This.Current_Mission.Trajectory.Nbody
    if i==1
        smint(i)=cspice_et2utc(This.Min_time(i),'C',0);
    else
        smint(i)=sprintf("%3.1f",round(This.Min_time(i)/60/60/24));
    end
end

%
% Now Construct listbox for Minimum Allowable Mission Times
%
mint=uicontrol('Style', 'listbox','Callback',@Min_Time_Set);
mint.FontSize=12;
mint.String=smint;
mint.Position=[580 230 120 290];

%
% Prepare String array for Maximum Allowable Mission Times for Optimization
%
smaxt=strings(This.Current_Mission.Trajectory.Nbody);

for i=1:This.Current_Mission.Trajectory.Nbody
    if i==1
        smaxt(i)=cspice_et2utc(This.Max_time(i),'C',0);
    else
        smaxt(i)=sprintf("%3.1f",round(This.Max_time(i)/60/60/24));
    end
end

%
% Now Construct listbox for Maximum Allowable Mission Times
%
maxt=uicontrol('Style', 'listbox','Callback',@Max_Time_Set);
maxt.FontSize=12;
maxt.String=smaxt;
maxt.Position=[720 230 120 290];

%
% Prepare String array for Initial Guess of Mission Times for Optimization
%
snomt=strings(This.Current_Mission.Trajectory.Nbody);

for i=1:This.Current_Mission.Trajectory.Nbody
    if i==1
        snomt(i)=cspice_et2utc(This.Current_Mission.Mission_Times(i),'C',0);
    else
        snomt(i)=sprintf("%3.1f",round(This.Current_Mission.Mission_Times(i)/60/60/24));
    end
end

%
% Construct Listbox for Initial Guess of Mission Times for Optimization
%
nomt=uicontrol('Style', 'listbox','Callback',@Nom_Time_Set);
nomt.FontSize=12;
nomt.String=snomt;
nomt.Position=[860 230 120 290];    

%
% Prepare String array for Minimum PERIAPSIS ALTITUDES for Optimization
% NOTE: IF The Body is an INTERMEDIATE POINT THEN THIS VALUE Should equal
% the distance of the POINT from the Sun in AU
%

sper=strings(This.Current_Mission.Trajectory.Nbody,2);

for i=1:This.Current_Mission.Trajectory.Nbody
    if This.Current_Mission.Trajectory.Body_Set(i).Fixed_Point>0
        sper(i,1)=sprintf("%12.6f",This.Min_Per(i)/This.AU);
        if (This.Max_dV(i)>=1e50)
            sper(i,2)="MAX";
        else
            sper(i,2)=sprintf("%12.6f",This.Max_dV(i)/1000);
        end
    else
        sper(i,1)=sprintf("%12.6f",This.Min_Per(i)/1000);
        if (This.Max_dV(i)>=1e50)
            sper(i,2)="MAX";
        else
            sper(i,2)=sprintf("%12.6f",This.Max_dV(i)/1000);
        end
    end
end

%
% Construct ENCOUNTER CONSTRAINT DETAILS Listbox
%
per=uicontrol('Style', 'listbox','Callback',@Enc_Con_Set);
per.FontSize=12;
per.String=sper(:,1);
per.Position=[1030 230 120 290]; 

%
% Construct ENCOUNTER CONSTRAINT DETAILS Listbox
%
per2=uicontrol('Style', 'listbox','Callback',@Enc_Con_Set);
per2.FontSize=12;
per2.String=sper(:,2);
per2.Position=[1200 230 120 290]; 
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Set_Mission wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Set_Mission_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%
% This function is executed if Planet name ListBox is Selected
%
function Planet_Name_Set(source,event)

global This;
global splan;
global plan;


val= source.Value;
info = source.String;

%
% Return if selection is not in required range
%

if (val > This.Current_Mission.Trajectory.Nbody)
    return;
end
 
    planinput=inputdlg("Enter Name of Body",sprintf("Planet %d",val),1,info(val),'on');
    if isempty(planinput)
        planinput = splan(val);
    end

splan(val)=planinput;
plan.String=splan;

for i=1:This.Current_Mission.Trajectory.Nbody
    This.Current_Mission.Trajectory.Body_Set(i).name=splan(i);
end

%
% This function is executed if Minimum Time Listbox is Selected
%
function Min_Time_Set(source,event)

global This;
global smint;
global mint;

val= source.Value;
info = source.String; 
%
% Return if selection is not in required range
%

if (val > This.Current_Mission.Trajectory.Nbody)
    return;
end

if val==1
    timeinput=inputdlg("Enter Launch Window Opening Date",sprintf("Planet %d",val),1,info(val));
else
    timeinput=inputdlg("Enter Minimum Cruise Time in Days (min = 1 day)",sprintf("Planet %d",val),1,info(val));
end

if isempty(timeinput)
    timeinput=smint(val);
end

for i=1:This.Current_Mission.Trajectory.Nbody
    if i == val
        if i==1
            try 
                dummy= cspice_str2et(char(timeinput));
            catch
                timeinput=smint(1);
            end
        else
            if isnan(str2double(timeinput))
                timeinput=smint(i);
            end
        end
        smint(i) = timeinput;
    end
    
    if i==1
        This.Min_time(i) = cspice_str2et(char(smint(i)));
    else
        if str2double(smint(i))<1.0
            smint(i)="1.0";
        end
        This.Min_time(i) = 24*60*60*str2double(smint(i));
    end
end
mint.String=smint;

%
% This function is executed if Maximum Time Listbox is Selected
%
function Max_Time_Set(source,event)

global This;
global smaxt;
global maxt;

val= source.Value;
info = source.String; 

%
% Return if selection is not in required range
%

if (val > This.Current_Mission.Trajectory.Nbody)
    return;
end

if val==1
    timeinput=inputdlg("Enter Launch Window Closing Date",sprintf("Planet %d",val),1,info(val));
else
    timeinput=inputdlg("Enter Maximum Cruise Time in Days (min = 1 day)",sprintf("Planet %d",val),1,info(val));
end

if isempty(timeinput)
    timeinput=smaxt(val);
end

for i=1:This.Current_Mission.Trajectory.Nbody
    if i == val
        if i==1
            try 
                dummy= cspice_str2et(char(timeinput));
            catch
                timeinput=smaxt(1);
            end
        else
            if isnan(str2double(timeinput))
                timeinput=smaxt(i);
            end
        end
        smaxt(i) = timeinput;
    end
    
    if i==1
        This.Max_time(i) = cspice_str2et(char(smaxt(i)));
    else
        if str2double(smaxt(i))<1.0
            smaxt(i)="1.0";
        end
        This.Max_time(i) = 24*60*60*str2double(smaxt(i));
    end
end
maxt.String=smaxt;

%
% This function is executed if Initial Guess Time Listbox is Selected
%
function Nom_Time_Set(source,event)

global This;
global snomt;
global nomt;

val= source.Value;
info = source.String; 

%
% Return if selection is not in required range
%

if (val > This.Current_Mission.Trajectory.Nbody)
    return;
end

if val==1
    timeinput=inputdlg("Enter Nominal Launch Date",sprintf("Planet %d",val),1,info(val));
else
    timeinput=inputdlg("Enter Nominal Cruise Time in Days (min = 1 day)",sprintf("Planet %d",val),1,info(val));
end

if isempty(timeinput)
    timeinput=snomt(val);
end




for i=1:This.Current_Mission.Trajectory.Nbody
    if i == val
        if i==1
            try 
                dummy= cspice_str2et(char(timeinput));
            catch
                timeinput=snomt(1);
            end
        else
            if isnan(str2double(timeinput))
                timeinput=snomt(i);
            end
        end
        snomt(i) = timeinput;
    end
    
    if i==1
        This.Current_Mission.Mission_Times(i) = cspice_str2et(char(snomt(i)));
    else
        if str2double(snomt(i))<1.0
            snomt(i)="1.0";
        end
        This.Current_Mission.Mission_Times(i) = 24*60*60*str2double(snomt(i));
    end
end

nomt.String=snomt;

%
% This function is executed if Minimum Perihelia Listbox is Selected
%
function Perihelia_Set(source,~)

global This;
global sperhel;
global perh;

val= source.Value;
info= source.String;
%
% Return if selection is not in required range
%

if (val > This.Current_Mission.Trajectory.Nbody)
    return;
end
    perhelinput=inputdlg("Enter Min. Perihelion in AU, set to zero if undefined",sprintf("Transfer %d to %d",val-1,val),1,info(val));
    
    if isempty(perhelinput)
        perhelinput=sperhel(val);
    end

    sperhel(val)=perhelinput;
    
    This.Perihelia_flag = 0;

for i=1:This.Current_Mission.Trajectory.Nbody
    tempstr= sperhel(i);
    if contains(sperhel(i),"P")
        This.Perihelia_flag = This.Perihelia_flag + 2^(i-1);
        tempstr=erase(sperhel(i),"P");        
    elseif contains(sperhel(i),"p")
        This.Perihelia_flag = This.Perihelia_flag + 2^(i-1);
        tempstr=erase(sperhel(i),"p");
    end
    
    if (~isnan(str2double(tempstr))) 
        This.Perihelia(i)=This.AU*str2double(tempstr);
    else
        sperhel(i) = perh.String(i);
    end
    
end
perh.String=sperhel;

%
% This function is executed if ENCOUNTER CONSTRAINT Listbox is Selected
%
function Enc_Con_Set(source,~)

global This;
global sper;
global per;
global per2;

val= source.Value;
info = source.String; 

%
% Return if selection is not in required range
%

if (val > This.Current_Mission.Trajectory.Nbody)
    return;
end

for i=1:This.Current_Mission.Trajectory.Nbody
    sper(i,1)=per.String{i};
    sper(i,2)=per2.String{i};
end
if This.Current_Mission.Trajectory.Body_Set(val).Fixed_Point==0
    
    answer = {'',''};
    prompt = {'Enter Minimum Periapsis Altitude in km:','Enter Constraint on DeltaV in km/sec:'};
    dlg_title = sprintf("Body %d",val);
    num_lines = 1;
    defaultans{1}=sprintf('%12.6f',This.Min_Per(val)/1000);
    
    if (This.Max_dV(val)>=1e50)
        defaultans{2}='MAX';
    else
        defaultans{2}=sprintf('%12.6f',This.Max_dV(val)/1000);
    end
    
    while (isnan(str2double(answer{1}))||isnan(str2double(answer{2})))
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans,'on');
        if (isempty(answer))
            break;
        end
        if (answer{2}=="MAX")
            answer{2}="1e50";
        end
    end
    
    if(~isempty(answer))
        This.Max_dV(val)=str2double(answer{2})*1000;
        This.Min_Per(val)=str2double(answer{1})*1000;
    end
    % See if Initial Periapsis value is set
    
    if (This.Min_Per(1)>0)
        This.Current_Mission.home_periapsis=This.Min_Per(1)+This.Current_Mission.Trajectory.Body_Set(1).radius;
    else
        This.Current_Mission.home_periapsis=0.0;
    end
    
    % See if Target Periapsis value set
    if (This.Min_Per(This.Current_Mission.Trajectory.Nbody)>0)
        This.Current_Mission.target_periapsis=This.Min_Per(This.Current_Mission.Trajectory.Nbody)+This.Current_Mission.Trajectory.Body_Set(This.Current_Mission.Trajectory.Nbody).radius;
    else
        This.Current_Mission.target_periapsis=0.0;
    end
    
    sper(val,1)=sprintf("%12.6f",This.Min_Per(val)/1000);
    if (This.Max_dV(val)>=1e50)
        sper(val,2)="MAX";
    else
        sper(val,2)=sprintf("%12.6f",This.Max_dV(val)/1000);
    end
else
    
    answer = {'',''};
    prompt = {'ENTER DISTANCE OF INTERMEDIATE POINT FROM CENTRE OF ECLIPTIC IN AU:','Enter Constraint on DeltaV in km/sec:'};
    dlg_title = sprintf("Encounter Constraint Details for Body %d",val);
    num_lines = 1;
    defaultans{1}=sprintf('%12.6f',This.Min_Per(val)/This.AU);
    
    if (This.Max_dV(val)>=1e50)
        defaultans{2}='MAX';
    else
        defaultans{2}=sprintf('%12.6f',This.Max_dV(val)/1000);
    end
    
    while (isnan(str2double(answer{1}))||isnan(str2double(answer{2})))
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        if (isempty(answer))
            break;
        end
        if (answer{2}=="MAX")
            answer{2}="1e50";
        end
        
    end
    if (~isempty(answer))
        This.Max_dV(val)=str2double(answer{2})*1000;
        This.Min_Per(val)=str2double(answer{1})*This.AU;
    end
    sper(val,1)=sprintf("%12.6f",This.Min_Per(val)/This.AU);
    if (This.Max_dV(val)>=1e50)
        sper(val,2)="MAX";
    else
        sper(val,2)=sprintf("%12.6f",This.Max_dV(val)/1000);
    end

end


per.String=sper(:,1);
per2.String=sper(:,2);
%
% This function is executed if Maximum Mission Duration edit box is
% selected
%
function Max_Dur_Set(source,event)

global This;
global durstring;

val= source.Value;
info = source.String;

if (info=="MAX")
    This.Max_Duration=This.MAX_DURATION;
    durstring = "MAX";
else
    if (isnan(str2double(info)))
        info=sprintf("%5.2f", This.Max_Duration/365/24/60/60);
    else
        This.Max_Duration = 365*24*60*60*str2double(info);
    end
    if(This.Max_Duration >= This.MAX_DURATION)
        This.Max_Duration = This.MAX_DURATION;
        durstring="MAX";
    else
        durstring = info;
    end
end
        

%
% This function is executed if Maximum RunTime for Optimization is Selected
%
function Run_Time_Set(source,event)

global This;

val= source.Value;
info = source.String; 

if (isnan(str2double(info)))
    info = sprintf("%5.2f", This.Run_Time/60 );
end
This.Run_Time = 60*str2double(info);

%
% This function is executed if the Rendezvous with target checkbox is selected
%
function Destination_Arrival_Set(source,event)

global This;
val= source.Value;
info = source.String;

This.Current_Mission.FlybyRendez = val;


%
% This function is executed if the Prograde checkbox is selected
%
function Prograde_Set(source,event)

global This;
val= source.Value;
info = source.String;

This.Current_Mission.wayflag = val;

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
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
