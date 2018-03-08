addpath Bodies;
addpath Transfers;
addpath Transfers/Gravity_Assist;
addpath Transfers/Gravity_Assist/With_Planetary_Encounter_Dynamics;

% Astronomical Unit in Metres
AU = 149597870700;


% Create Project to Work on in OITS

global This;

This = Project;

This.Body_Select=Body(This.Max_NBody);

This.NBody_List=0;

% Initialize NLopt Optimization Software

addpath('NLOPT\matlab');
addpath('NLOPT\');

mex -output NLOPT\nlopt_optimize -LNLOPT -lnlopt-0 -INLOPT\ NLOPT\matlab\nlopt_optimize.c
  
% Use Default Name and BSP file (Binary SPK file from NASA)

This.name = 'Test';

This.BSP = 'SPICE\de430.bsp';

% Initialize Pointer to List of Bodies

global Body_pointer;

Body_pointer=1;
This.Body_Number=Body_pointer;

This = This.Initialize_SPICE;

This = This.Get_SPICE_List(This.BSP);
This = This.Add_Intermediate_Point;
This = This.Add_Fixed_Point;
This = This.Get_SPICE_List('SPICE\1000012.bsp');
This = This.Get_SPICE_List('SPICE\lutetia.bsp');
This = This.Get_SPICE_List('SPICE\steins.bsp');
This = This.Get_SPICE_List('SPICE\extrasolar.bsp');
This = This.Get_SPICE_List('SPICE\101955.bsp');
This = This.Merge_Body_Data;

% Initialize Selected Body List to First of the Bodies Available
This.Body_Select(1)=This.Body_List(1);

f=figure(OITS);

% txt=uicontrol(f,'Style','text','Position', [150 10 250 40],'String',This.name);
uiwait(f);


This

%f.Visible='on';
