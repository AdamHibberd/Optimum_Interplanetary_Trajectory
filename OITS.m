function varargout = OITS(varargin)
% OITS MATLAB code for OITS.fig
%      OITS, by itself, creates a new OITS or raises the existing
%      singleton*.
%
%      H = OITS returns the handle to a new OITS or the handle to
%      the existing singleton*.
%
%      OITS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OITS.M with the given input arguments.
%
%      OITS('Property','Value',...) creates a new OITS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OITS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OITS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OITS

% Last Modified by GUIDE v2.5 20-Oct-2017 09:15:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OITS_OpeningFcn, ...
                   'gui_OutputFcn',  @OITS_OutputFcn, ...
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


% --- Executes just before OITS is made visible.
function OITS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OITS (see VARARGIN)

% Choose default command line output for OITS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OITS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OITS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
%
% This Function Executes on Selection of OPEN MISSION in Main Menu
%
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global This;

% Construct default filename from Project name (Use .mat matlab filetype)
default_file = strcat(This.name,'.mat') ; 

% Get File
filename = uigetfile(default_file);

if (filename == 0 )
    
else
    
    % Load Project Object from File into global Project 'This'
    
    load(filename, 'This');
    
    % Old Style Files Convert
    if(size(This.Max_dV)<=1)
        This.Max_dV(1:This.Current_Mission.Trajectory.Nbody)=1e50;
        if(This.Min_Per(1)>0)
            This.Max_dV(1)=sqrt(This.Min_Per(1))*1000;
            This.Min_Per(1)=0.0;
        end
    end
end


% --- Executes on button press in pushbutton2.
%
% This Function Executes on selection of NEW MISSION in main Menu
%
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deal with selection of planets for This Mission

global This;            % This Mission
global Body_pointer;    % Pointer to Current Location in List of Bodies
global Itemstemp;       % String List of Selected Bodies
global txt;             % Tag for text box
global lst;             % Tag for list box


Itemstemp=cell(This.Max_NBody);

%
% Make Figure to allow selection of available Solar Systen Objects As Read
% in from SPICE KERNELS
%
%

figureselect=figure;
figureselect.ToolBar='none';
figureselect.MenuBar='none';
figureselect.Name='Select From Available Solar System Objects';
pos=figureselect.Position;
figureselect.Position= [300   378   900   420];


%
% Initialize array of Names of Available Solar System Objects to be
% displayed in the popupmenu
%
%

for i=1:This.NBody_List
    Items(i)=cellstr(This.Body_List(i).name);
end

%
% Initialize List of Selected Bodies
%
%

for i=1:This.Body_Number
    Itemstemp(i)=cellstr(This.Body_Select(i).name);
end


%
% Construct popupmenu in the current figure on LHS. These are all the
% available SPICE SSO's from which to choose
%
%

pop=uicontrol(figureselect,'Style', 'popupmenu', 'String', Items, 'Callback',@SelectBody);
pop.Position = [125 120 300 200];
txt=uicontrol(figureselect,'Style','text');
txt.Position=[110  320 300 50];
txt.FontSize=14;
txt.String=sprintf("Select Body Number %d",Body_pointer);

%
% Construct listbox in the current figure on RHS of currently selected SSO's from popupmenu
%
%
lst=uicontrol(figureselect,'Style', 'listbox', 'String', Itemstemp,'Callback', @GoBack);
lst.Position=[500 120 300 200];

%
% Construct pushbutton above listbox to allow exit of current Figure 
%
%
but=uicontrol(figureselect,'Style', 'pushbutton', 'String', 'Quit Selection','Callback',@QuitSelect);
but.Position=[500 340 180 40];
but.FontSize=14;

%
% Construct pushbutton above listbox to allow clearance of currently selected list of SSO's
%
%

cle=uicontrol(figureselect,'Style', 'pushbutton', 'String', 'Clear Selection','Callback',@ClearSelect);
cle.Position=[700 340 180 40];
cle.FontSize=14;


% --- Executes on button press in pushbutton3.
%
% This Function Executes on Selection of SAVE MISSION in Main Menu
%
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global This;

%
% Save the current Project object 'This' into user specified file
%
%
uisave('This',This.name);

    

% --- Executes on button press in pushbutton4.
%
% This Function Executes on Selection of SET OPTIMIZER DETAILS in Main Menu
% It Opens a new figure Called SET_MISSION which allows more detailed
% Project Settings
%

function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global This;

%
% Only Open Set_Mission if at least two SSO's have been selected
%
%
if (This.Body_Number>1)     
        This.factor =0.5;
        
        %
        % Check to see if Set_Mission has already been selected
        %
        if(isempty(This.Current_Mission))   
        
                %
                % Check to see if The Optimizer has been executed on this mission yet
                %
                if(isempty(This.Solution))  
                    
                    %
                    % If this is a new mission construct List(array) of chosen Body Objects
                    %
                    This.Body_Chosen = Body(This.Body_Number);
                    
                    %
                    % Inititalise List of chosen bodies
                    %
                    for i=1:This.Body_Number                    
                        This.Body_Chosen(i) = This.Body_Select(i);
                        
                        %
                        % Inialise array of times
                        %
                        if (i==1)
                            time(i)=cspice_str2et(char(datetime('today')));     
                        else
                            time(i)=time(i-1)+365*24*60*60;
                        end
                    end
                    
                    %
                    % If this is a new mission, Construct the Mission based
                    % on data obtaines above
                    %
                    This.Current_Mission = Mission( This.Body_Chosen,time,This.Min_Spice_Select,This.Max_Spice_Select);
                else
                    
                    %
                    % Set Current Mission to Solution of Optimizer Run
                    % (This code is precautionary)
                    %
                    This.Current_Mission = This.Solution;
                end
                            
        end
        
        %
        % Create The Set_Mission Figure based on details derived above or on
        % The Current_Mission details as last specified
        %
        
        g = figure(Set_Mission);
        
        uiwait(g); 
        
        l=gcf;
        l.Visible ='off';
        
        % Ensure That Minimum times are all less than maxima and
        % Ensure Initial Guess Lies in within min and max times

        
        while max((This.Min_time>This.Current_Mission.Mission_Times)|(This.Max_time<This.Current_Mission.Mission_Times))>0
                  h=warndlg('ERROR FOUND IN TIME BOUNDS OR INITIAL GUESS');
                  uiwait(h);
                  k = figure(Set_Mission);
                  uiwait(k);         
        end

        l.Visible='on';
        
        % Initialise Mission
  %      if (isempty(This.Solution))
  %          
  %          FlybyRendez=This.Current_Mission.FlybyRendez;
  %          wayflag=This.Current_Mission.wayflag;
  %          This.Current_Mission = This.Current_Mission.Set_Absolute_Times( This.Current_Mission.Mission_Times );
  %          This.Current_Mission = Mission( This.Body_Chosen,This.Current_Mission.Absolute_Times,This.Min_Spice_Select,This.Max_Spice_Select);
  %          This.Current_Mission.FlybyRendez=FlybyRendez;
  %          This.Current_Mission.wayflag=wayflag;
%
 %       end
end

% --- Executes on button press in pushbutton5.
%
% This Function Executes on Selection of VIEW RESULTS in Main Menu
%
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global This;

%
% Firstly Display Basic Info in numeric form about the encounter with each SSO in turn
%
This = This.View_Info(2);

%
% Now Display Interplanetary Trajectory in 2D & 3D Plot form as well as orbits
% of SSO's. Also Display Speed Against Radial Distance Plot of Trajectory
%
This = This.View_Results( 600, 2);

%
% Now Display Orbital Information in numeric form
%
This = This.View_Orbit_Info( 2 );

%
% Now Display Plots of Distance Against Time for Each SSO Encounter as well
% As 3D Trajectory for each SSO Encounter
%
if (This.Body_Number>2)
    This = This.View_Planetary_Encounters( 600, 2 );
end

 %This = This.View_DeltaV_Vs_Time(40,2,365*24*60*60);


% --- Executes on button press in pushbutton6.
%
% This Function Executes on Selection of RUN OPTIMIZATION in Main Menu
%
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global This;

%
% Only Optimize if At least 2 Objects have been Selected
%
if (This.Body_Number>1)
    
    %
    % Initialize The Perihelia Constraints if They have been Chosen
    %
    This.NPerihelia = 0;
    This.Nconstraints = This.Current_Mission.Trajectory.Nbody-2;
    for i=2:This.Current_Mission.Trajectory.Nbody
        if (This.Perihelia(i)>0.0)
            This.NPerihelia = This.NPerihelia+1;
            This.Per_Pointer(This.NPerihelia)= i;
        end
    end
    
    %
    % Run Optimizer On Current Project / Mission Combination
    %
    This = This.Optimize_Mission(4);

end

%
% This Function Executes on Selection of a Body in the popupmenu on the 
% LHS of NEW MISSION Figure
%
function SelectBody(source,event)
    
    global This;
    global Body_pointer;    % Initialized to 
    global Itemstemp;
    global txt;
    global lst;
   
    val= source.Value;      % This is the row number of the popup selected
    info = source.String;   % This is not used
    
    %
    % New Mission Needs to be generated, any old missions are deleted
    %
    This.Current_Mission=[];
    This.Solution=[];
    
    %
    % Add the selected Body on the LHS to the Selected Body List with the
    % respective minimum and maximum SPICE times
    %
    This.Body_Select(Body_pointer)=This.Body_List(val);
    This.Min_Spice_Select(Body_pointer) = This.Min_Spice_Time(val);
    This.Max_Spice_Select(Body_pointer) = This.Max_Spice_Time(val);
    
    %
    % Number of Selected Bodies (Already Initialized to 1 in run_OITS)
    %
    This.Body_Number = max(Body_pointer,This.Body_Number);
    
    Body_pointer = Body_pointer+1;
    txt.String=sprintf("Select Body Number %d",Body_pointer);
    
    %
    % Display String for Selected Body List on RHS of Figure
    %
    for i=1:This.Body_Number
        Itemstemp(i)=cellstr(This.Body_Select(i).name);
    end
    lst.String=Itemstemp;
    
%
% This function is executed if User Selects one of the SSO's that he has
% already selected on the RHS of the NEW MISSION figure
%
function GoBack(source,event)

    global Body_pointer;
    global txt;
    val= source.Value;
    info = source.String;
    
    Body_pointer=val;
    txt.String=sprintf("Select Body Number %d",Body_pointer);

%
% This function is executed if the Quit Option Is selected in the NEW
% MISSION figure
%
function QuitSelect(source,event)

global This;

%
% Only Initialise Current Mission if more than one SSO is Selected
%
if (This.Body_Number>1)
    
    %
    % Construct the Body list of Chosen Objects 
    %
    This.Body_Chosen = Body(This.Body_Number);
     
    %
    % Construct and zero the List of Time Ranges needed for the Optimizer,
    % THESE MUST LATER BE SPECIFIED BY USER
    %
    This.Min_time=zeros(1,This.Body_Number);
    This.Max_time=zeros(1,This.Body_Number);
    
    %
    % Construct and zero:
    % 1) List of Minimum Periapsis Distances
    % 2) List of Minimum Perihelia Disances
    % THESE MUST LATER BE SPECIFIED BY USER
    %
    This.Min_Per=zeros(1,This.Body_Number);
    This.Max_dV=zeros(1,This.Body_Number);
    This.Perihelia=zeros(1,This.Body_Number);

    % The following must be specified by user.
    for i=1:This.Body_Number
        This.Min_Per(i)=0.0;
        This.Max_dV(i)=1e50;
    end
    
    %
    % Initialize the Body List of Chosen Objects
    %
    for i=1:This.Body_Number
        This.Body_Chosen(i) = This.Body_Select(i);
        if (i==1)
            time(i)=cspice_str2et(char(datetime('today')));
        else
            time(i)=time(i-1)+2*365*24*60*60;
        end
    %
    % Initialize the Fixed Points
    %
        if This.Body_Chosen(i).Fixed_Point>0
            This.Body_Chosen(i).ephem0.r = [This.AU*cos(i*pi/100) This.AU*sin(i*pi/100) 0];
            This.Body_Chosen(i).ephem0.v = [ 0 0 0 ];
            This.Body_Chosen(i).ephem0.t = time(i);
            This.Min_Per(i)=norm(This.Body_Chosen(i).ephem0.r);
        end
    end
    
    %
    % Now Construct The Current Mission Object and Initialize as Required
    %
    This.Current_Mission = Mission( This.Body_Chosen,time,This.Min_Spice_Select,This.Max_Spice_Select);

    for i=1:This.Body_Number
        This.Min_time(i)=This.Current_Mission.Mission_Times(i)-60*60*24*28;
        This.Max_time(i)=This.Current_Mission.Mission_Times(i)+60*60*24*28;
    end
end

close;

%
% This function is executed if the Clear Selection Option Is selected in the NEW
% MISSION figure
%

function ClearSelect(source,event)

global This;
global Body_pointer;
global txt;
global lst;
global Itemstemp;

This.Body_Select=Body;

for i=1:This.Body_Number
    Itemstemp(i)={""};
end

Body_pointer=1;

This.Body_Number=1;  
   
    txt.String=sprintf("Select Body Number %d",Body_pointer);
    lst.String=Itemstemp;
 
   
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    close;


% --- Executes on button press in pushbutton8.
%
% This Function Executes on Selection of ANIMATE RESULTS in Main Menu
%
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global This;

title=inputdlg("Enter Title for Animation","");

% titleanim=input('Enter Title of Animation:','s');
% This= This.Animate_Results(int64(5000/This.Current_Mission.Trajectory.Nbody), 2, titleanim);
 if (~isempty(title))
     titleanim=char(title);
     w1=msgbox('Starting Animation Please Wait - You Will be Informed When Animation is Complete','');
     
     This= This.Animate_Results(int64(5000/This.Current_Mission.Trajectory.Nbody), 2, titleanim);
     w2=helpdlg('Animation Complete and Stored in File TrajVideo.mp4','');
     uiwait(w2);
     
 end
