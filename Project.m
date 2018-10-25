%# Project is a class for Optimizing Interplanetary Trajectory Missions as specifed by user
classdef Project
%#
%# PROPERTIES:
%#
%#        
%#    name;               : Name of Project
%#    BSP;                : Filename for Binary SPK file from JPL
%#    Body_List;          : List of Bodies Available to Choose from
%#    NBody_List;         : Number of Bodies Available to Choose from
%#    Max_NBody = 20;     : Maximum Number of Bodies allowed for Optimizer
%#    Min_Spice_Time;     : Lower Limit of Spice Kernel Range
%#    Max_Spice_Time;     : Upper Limit of Spice Kernel Range
%#    Min_Spice_Select;   : Selected Values of Min_Spice_Time
%#    Max_Spice_Select;   : Selected Values of Max_Spice_Time
%#    Body_Select;        : List of Selected Bodies
%#    Body_Chosen;        : List of Chosen Bodies for Optimization
%#    Body_Number;        : Number of Selected Body
%#    Current_Mission;    : Current Mission Selected bu User
%#    Global_Solution;    : Result of Running Global Optimizer
%#    Local_Solution;     : Result of Running Local Optimzer
%#    Solution;           : Result of Last run of Optimizer
%#    Run_Time = 60;      : Maximum Optimizer Run Time (Default 60secs)
%#    Nconstraints;       : Number of Periapsis Constraints
%#    NPerihelia;         : Number of Perihelia Constraints
%#    Per_Pointer;        : Pointer to Array of Perihelia Values
%#    Constr_Tol;         : Array of Constraint Tolerances
%#    MAX_DURATION=1e50;  : Upper Limit for Max_Duration
%#    Max_Duration=1e50;  : Maximum Duration Of Entire Mission (optional)
%#    factor = 0.5;       : Factor used for Local Optimizer
%#    AU = 149597870700;  : Astronomical Unit in Metres
%#    Min_time;           : Array of Minimum Times for Optimization
%#    Max_time;           : Array of Maximum Times for Optimization
%#    Min_Per;            : Array of Minimum Periapsis for each Body
%#    Max_dV;             : Array of Maximum Allowable DeltaV's for each Body
%#    Perihelia;          : Array of Minimum Perihelia for each Transfer
%#    Perihelia_flag=0;   : Flag to indicate if Perihelia Constraints are Orbital Parameters or achieved
%#
%# METHODS:
%#
%# Initialize_SPICE         :   Initializes SPICE Toolkit and Opens Leap Second File naif0012.tls    
%# Get_SPICE_List           :   Opens and extracts BINARY SPICE KERNEL Files .BSP
%# Merge_Body_Data          :   Merges Data on the Planets with Data extracted from Get_SPICE_List  
%# Add_Intermediate_Point   :   Adds an INTERMEDIATE POINT to the Possible Bodies to Select From
%# Add_Fixed_Point          :   Adds a FIXED POINT to the Possible Bodies to Select From
%# Optimize_Mission         :   Optimizes Mission provided by Current_Mission -> Solution goes to Solution
%# Compute_DeltaV_NLopt     :   Calculates DeltaV as determined by times input by Optimize_Mission
%# Update_Traj              :   Updates current Trajectory and calculates all the constraints based on the times and thetas and thi's
%# Per_NLopt                :   Calculates Periapsis Constraints for Optimize_Mission
%# dV_NLopt                 :   Calculates DeltaV Constraints for Optimize_Mission
%# Perhel                   :   Calculates Perihelion Constraints for Optimize_Mission
%# Overall_Duration         :   Calculates Overall Mission Duration Constraints for Optimize_Mission
%# View_Results             :   Displays Interplanetary Trajectory in 2D & 3D Plot form as well as orbits
%#                          :   of SSO's. Also Displays Speed Against Radial Distance Plot of Trajectory
%# View_Info                :   Display Basic Info in numeric form about the encounter with each SSO in turn
%# View_Planetary_Encounters:   Displays Plots of Distance Against Time for Each SSO Encounter as well
%#                          :   As 3D Trajectory for each SSO Encounter
%# View_DeltaV_Vs_Time      :   Displays Plots of DeltaV Against Time
%# View_Orbit_Info          :   Displays Orbital Information in numeric form
%# Animate_Results          :   Animates Interplanetary Trajectory
%#

properties
        
    name;               % Name of Project
    BSP;                % Filename for Binary SPK file from JPL
    Body_List;          % List of Bodies Available to Choose from
    NBody_List;         % Number of Bodies Available to Choose from
    Max_NBody = 20;     % Maximum Number of Bodies allowed for Optimizer
    Min_Spice_Time;     % Lower Limit of Spice Kernel Range
    Max_Spice_Time;     % Upper Limit of Spice Kernel Range
    Min_Spice_Select;   % Selected Values of Min_Spice_Time
    Max_Spice_Select;   % Selected Values of Max_Spice_Time
    Body_Select;        % List of Selected Bodies
    Body_Chosen;        % List of Chosen Bodies for Optimization
    Body_Number;        % Number of Selected Body
    Current_Mission;    % Current Mission Selected bu User
    Global_Solution;    % Result of Running Global Optimizer
    Local_Solution;     % Result of Running Local Optimzer
    Solution;           % Result of Last run of Optimizer
    Run_Time = 60;      % Maximum Optimizer Run Time (Default 60secs)
    Nconstraints;       % Number of Periapsis Constraints
    NPerihelia;         % Number of Perihelia Constraints
    Per_Pointer;        % Pointer to Array of Perihelia Values
    Constr_Tol;         % Array of Constraint Tolerances
    MAX_DURATION=1e50;  % Upper Limit for Max_Duration
    Max_Duration=1e50;  % Maximum Duration Of Entire Mission (optional)
    factor = 0.5;       % Factor used for Local Optimizer
    AU = 149597870700;  % Astronomical Unit in Metres
    Min_time;           % Array of Minimum Times for Optimization
    Max_time;           % Array of Maximum Times for Optimization
    Min_Per;            % Array of Minimum Periapsis for each Body
    Max_dV;             % Array of Maximum Allowable DeltaV for each Body
    Perihelia;          % Array of Minimum Perihelia for each Transfer
    Perihelia_flag=0;   % Flag to indicate if Perihelia Constraints are Orbital Parameters or achieved

end
    
methods

     function obj = Initialize_SPICE(obj)
        
%# Initialize_SPICE         :   Initializes SPICE Toolkit and Opens Leap Second File naif0012.tls        
        % Initialise Various SPICE files 
        
        addpath('thirdparty\SPICE\');
        addpath('thirdparty\SPICE\mice\mice\src\mice');
        addpath('thirdparty\SPICE\mice\mice\lib');
        cspice_tkvrsn('toolkit');
        cspice_furnsh('thirdparty\SPICE\naif0012.tls');
        
    end
 
    function obj = Get_SPICE_List(obj, SPK)
%# Get_SPICE_List           :   Opens and extracts BINARY SPICE KERNEL Files .BSP

      MAXIV  = 1000;
      WINSIZ = 2 * MAXIV;
        
      % Initialise Kernel
      cspice_furnsh(SPK);
      %
      % Find the set of objects in the SPK file.
      %
      ids = cspice_spkobj( SPK, MAXIV );
        
      BODYLISTN = Body(numel(ids));
      

          
      %
      % We want to display the coverage for each object. Loop over
      % the contents of the ID code set, find the coverage for
      % each item in the set, and display the coverage.
      %
      for i=1:numel(ids)

         %
         % Extract the coverage data for object 'ids(i)'.
         %
         cover     = cspice_spkcov( SPK, ids(i), WINSIZ );
         %
         % Display a simple banner.
         %
         fprintf( '========================================\n')
         fprintf( 'Coverage for object %d\n', ids(i) )

         %
         % Convert the endpoints to TDB calendar format time strings
         % and display them. Pass the endpoints in an array,
         % so cspice_timout returns an array of time strings.
         %
         % Recall a vectorized input has dimension 1xM so transpose
         % the 'cover' slice.
         %
         timstr = cspice_timout( cover(1:2)', ...
                            'YYYY MON DD HR:MN:SC.### (TDB) ::TDB' );
         fprintf('   Start: %s\n'  , timstr(1,:) )
         fprintf('    Stop: %s\n\n', timstr(2,:) )

         BODYLISTN(i).ID=sprintf('%d', ids(i));
         BODYLISTN(i).name = cspice_bodc2n(ids(i));
         if(BODYLISTN(i).name=="")
             BODYLISTN(i).name=BODYLISTN(i).ID;
         end
         BODYLISTN(i).radius=0;
         obj.Min_Spice_Time(i+obj.NBody_List) = cspice_str2et(timstr(1,:));
         obj.Max_Spice_Time(i+obj.NBody_List) = cspice_str2et(timstr(2,:));
      end
      
      if (obj.NBody_List==0)
          
          obj.Body_List = BODYLISTN;
          obj.NBody_List = numel(BODYLISTN);
          
      else
          
          obj.Body_List = cat (2, obj.Body_List, BODYLISTN);
          obj.NBody_List = numel(obj.Body_List);
      
      end
      
    end


    function obj = Merge_Body_Data( obj )
%# Merge_Body_Data          :   Merges Data on the Planets with Data extracted from Get_SPICE_List

        PLANET_NAMES = {"MERCURY", "VENUS", "EARTH", "MARS", "JUPITER", "SATURN", "URANUS", "NEPTUNE", "PLUTO"};
        PLANET_MU    = 6.67259e-11*[ 0.3302e24 4.869e24 5.9736e24 0.6419e24 1898.6e24 568.46e24 86.83e24 102.43e24 0.0125e24];
        PLANET_RADIUS= [2440e3 6052e3 6378135 3397e3 71492e3 60268e3 25559e3 24766e3 1137e3];

        for i=1:obj.NBody_List
            for j=1:9
                if(~isempty(strfind(obj.Body_List(i).name, PLANET_NAMES{j})))
                    obj.Body_List(i).radius=PLANET_RADIUS(j);
                    obj.Body_List(i).mu=PLANET_MU(j);
                    break;
                end
            end
            if (strcmp(obj.Body_List(i).ID,'3788040'))
                obj.Body_List(i).name="OUMUAMUA";   % Make sure OUMUAMUA has a name!
            end
        end
        
    end


    function obj= Add_Intermediate_Point( obj)
%# Add_Intermediate_Point   :   Adds an INTERMEDIATE POINT to be Optimized to the Possible Bodies to Select From

        Joker = Body;
        Joker.ID = 'INTERMEDIATE POINT';
        Joker.name = "INTERMEDIATE POINT";
        Joker.radius = 0;
        Joker.Fixed_Point = 1;
        Joker.ephem0.r = [ obj.AU 0 0 ];
        Joker.ephem0.v = [ 0 0 0 ];
        Joker.ephem0.t = 0;
        spice_min= -1e50;
        spice_max=  1e50;
        obj.Body_List = cat(2, obj.Body_List, Joker);
        obj.NBody_List = numel(obj.Body_List);
        obj.Min_Spice_Time(obj.NBody_List) = spice_min;
        obj.Max_Spice_Time(obj.NBody_List) = spice_max;
        
    end

     function obj= Add_Fixed_Point( obj)
%# Add_Fixed_Point   :   Adds a FIXED POINT to the Possible Bodies to Select From

        Joker2 = Body;
        Joker2.ID = 'FIXED POINT';
        Joker2.name = "FIXED POINT";
        Joker2.radius = 0;
        Joker2.Fixed_Point = 2;
        Joker2.ephem0.r = [ obj.AU 0 0  ];
        Joker2.ephem0.v = [ 0 0 0 ];
        Joker2.ephem0.t = 0;
        spice_min= -1e50;
        spice_max=  1e50;
        obj.Body_List = cat(2, obj.Body_List, Joker2);
        obj.NBody_List = numel(obj.Body_List);
        obj.Min_Spice_Time(obj.NBody_List) = spice_min;
        obj.Max_Spice_Time(obj.NBody_List) = spice_max;
        
    end
    function obj = Optimize_Mission(obj, ~)
%# Optimize_Mission         :   Optimizes Mission provided by Current_Mission -> Solution goes to Solution           
        
% Set up Inputs To Optimizer
      
    obj.Solution=obj.Current_Mission;
    
% Calculate number of Intermediate Points in Body_Set
    NIP =0;
    for i1=1:obj.Solution.Trajectory.Nbody
        if (obj.Solution.Trajectory.Body_Set(i1).Fixed_Point==1)
            NIP=NIP+1;
        end
    end
    
    tin=zeros(1,obj.Solution.Trajectory.Nbody+2*NIP);
    lb=zeros(1,obj.Solution.Trajectory.Nbody+2*NIP);
    ub=zeros(1,obj.Solution.Trajectory.Nbody+2*NIP);

    for i1=1:obj.Solution.Trajectory.Nbody
        lb(i1)=obj.Min_time(i1);
        ub(i1)=obj.Max_time(i1);
    end

    % Set up Inequality Constraint Functions
    funcstring= '@(x)[';


    % Firstly The Minimum allowable Periapsis
    if (obj.Nconstraints>0)
        for i1=1:obj.Nconstraints
            bb_output_type{i1}='PB';
            funcstring= strcat(funcstring, sprintf(' Per_NLopt(x,%d)',i1+1));
                 funcstring = strcat( funcstring ,' ;' );
            nle(i1)=-1;
            nlrhs(i1)=0;
        end
    end

    % Secondly The Maximum allowable Periapsis (This adds a scale to the
    % Minimum Perapsis and reduces the chance of overshoot of the minimum
    % periapsis.)
    
    if (obj.Nconstraints>0)
        for i1=(obj.Nconstraints+1):2*obj.Nconstraints
            bb_output_type{i1}='PB';
            funcstring= strcat(funcstring, sprintf(' Per_NLopt(x,%d)',i1+1));
            if i1<2*obj.Nconstraints
                 funcstring = strcat( funcstring ,' ;' );
            end
            nle(i1)=-1;
            nlrhs(i1)=0;
        end
    end

    % Thirdly The Minimum allowable Perihelia
    if(obj.NPerihelia>0)
        funcstring = strcat( funcstring ,' ;' );
        for i1=1:obj.NPerihelia
            funcstring = strcat(funcstring,sprintf(' Perhel(x,%d)',i1));
            if i1<obj.NPerihelia
                funcstring = strcat( funcstring ,' ;' );
            end
            nle(i1+2*obj.Nconstraints)=-1;
            nlrhs(i1+2*obj.Nconstraints)=0;
            bb_output_type{i1+2*obj.Nconstraints}='EB';
    % bb_output_type{i+2*obj.Nconstraints}='PB';
        end
    end

    Max_Dur_Flag = 0;
    % Check to see if limit on the Maximum Duration
    if (obj.Max_Duration < obj.MAX_DURATION )
        Max_Dur_Flag = 1;
        funcstring = strcat( funcstring ,' ; Overall_Duration(x)');
        nle(2*obj.Nconstraints+obj.NPerihelia+1)=-1;
        nlrhs(2*obj.Nconstraints+obj.NPerihelia+1)=0;
        bb_output_type{2*obj.Nconstraints+obj.NPerihelia+1}='EB';
    end

    % Finally check for Maximum DeltaV Constraint
    Number_DeltaV_Constraints=0;
    for i1=1:obj.Solution.Trajectory.Nbody
        if (obj.Max_dV(i1)<1e50)
            Number_DeltaV_Constraints = Number_DeltaV_Constraints + 1;
            funcstring = strcat(funcstring,sprintf(' ; dV_NLopt(x,%d)',i1));
            nle(2*obj.Nconstraints+obj.NPerihelia+Max_Dur_Flag+Number_DeltaV_Constraints)=-1;
            nlrhs(2*obj.Nconstraints+obj.NPerihelia+Max_Dur_Flag+Number_DeltaV_Constraints)=0;
            bb_output_type{2*obj.Nconstraints+obj.NPerihelia+Max_Dur_Flag+Number_DeltaV_Constraints}='PB';
        end
    end

    % Construct Function Handle From String.
    funcstring=strcat(funcstring, ' ]');
    nlcon = eval(funcstring);
    
    % Specify Problem type (Non-linear)
    
    optiSolver('NLP');
    
    % Initial Guess           
    
     % Check for Presence of Intermediate Point
     NIP=-1;
     for i1=1:obj.Solution.Trajectory.Nbody
         if obj.Solution.Trajectory.Body_Set(i1).Fixed_Point==1  % INTERMEDIATE POINT ?
             NIP=NIP+2;
             
             % Set-up initial values for the Ecliptic polar co-ordinates,
             % theta and phi for this INTERMEDIATE POINT
             
             R=norm(obj.Solution.Trajectory.Body_Set(i1).ephemt.r);
             theta= atan2(obj.Solution.Trajectory.Body_Set(i1).ephemt.r(2),obj.Solution.Trajectory.Body_Set(i1).ephemt.r(1));
             phi=asin(obj.Solution.Trajectory.Body_Set(i1).ephemt.r(3)/R);
             
             % Set-up initial value of tin
             
             if NIP==1
                tin = cat ( 2, obj.Solution.Mission_Times ,[ theta phi ] );
             else
                tin = cat ( 2, tin , [ theta phi ] );
             end

             % The lower and upper bounds on the theta and phi's
             
             lb(obj.Solution.Trajectory.Nbody+NIP)=-pi;
             lb(obj.Solution.Trajectory.Nbody+NIP+1)=-pi/2;
             ub(obj.Solution.Trajectory.Nbody+NIP)=pi;
             ub(obj.Solution.Trajectory.Nbody+NIP+1)=pi/2;                     
         end
     end
     
     % Make sure tin is initialised if no intermediate points!
    
     for i1=1:obj.Solution.Trajectory.Nbody
        tin(i1) = obj.Solution.Mission_Times(i1);
     end 

    % Initialise the values of the last updated variables used by the
    % Objective Function and the Constraint Functions
    
    tlast1 = tin - tin;
    DeltaVold=0;
    ceq=zeros(1,min(1,2*obj.Nconstraints));
    ceqold=zeros(1,min(1,2*obj.Nconstraints));
    dVeq=zeros(1,min(1,2*obj.Nconstraints));
    dVeqold=zeros(1,min(1,2*obj.Nconstraints));
    req=zeros(1,obj.Solution.Trajectory.Nbody-1);
    reqold=zeros(1,obj.Solution.Trajectory.Nbody-1); 

    % Initialise Trajectory

    DeltaV =  Compute_DeltaV_NLopt(tin);            
    
    % Initialise NOMAD Optimising Settings
    if (obj.Nconstraints>0||obj.NPerihelia>0)
        nopts = nomadset('bb_output_type',bb_output_type ,'vns_search',0.95,'max_eval',60000);
    elseif obj.Max_Duration < obj.MAX_DURATION
        nopts = nomadset('bb_output_type',bb_output_type,'vns_search',0.75);
    else
        nopts=nomadset('vns_search',0.75);
        nlrhs=[];
        nle=[];
    end
    opts=optiset('solver','nomad','display','iter','maxfeval',60000,'maxtime',obj.Run_Time,'solverOpts',nopts); 
    Opt = opti('fun',@Compute_DeltaV_NLopt,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'x0',tin,'options',opts);

    % Run Optimization
    
    [Optimum,DeltaV,~,~] = solve(Opt,tin) 
   
 %   tin
 %   Optimum

    tlast1 = tin - tin;

    DeltaV          

    tin = Optimum;
    
    %Set Up the Trajectory for storing
    
    for i1=1:obj.Solution.Trajectory.Nbody
       obj.Solution.Mission_Times(i1) = Optimum(i1);
    end
    [DeltaV,~] =  Compute_DeltaV_NLopt(tin)

    % Set Solution 
    obj.Current_Mission = obj.Solution;
    obj.Global_Solution = obj.Solution;
    
    return;

        function [DeltaV, gradient  ] = Compute_DeltaV_NLopt(tin)
%# Compute_DeltaV_NLopt     :   Calculates DeltaV as determined by times input by Optimize_Mission

           if ~isequal(tin,tlast1)
                [DeltaV,gradient] = Update_Traj(tin);
                DeltaVold=DeltaV;
           end

            DeltaV=DeltaVold;
            return;
        end
       
        function [DeltaV,gradient] = Update_Traj(tin)
%# Update_Traj     :   Updates current Trajectory and calculates all the constraints based on the times and thetas and thi's i.e. tin

            % Firstly calculate the Intermediate Points 
            NIP=-1;
            for j=1:obj.Solution.Trajectory.Nbody
                if obj.Solution.Trajectory.Body_Set(j).Fixed_Point==1
                    NIP=NIP+2;
                    coslong=cos(tin(obj.Solution.Trajectory.Nbody+NIP));
                    sinlong=sin(tin(obj.Solution.Trajectory.Nbody+NIP));
                    coslat=cos(tin(obj.Solution.Trajectory.Nbody+NIP+1));
                    sinlat=sin(tin(obj.Solution.Trajectory.Nbody+NIP+1));
                    obj.Solution.Trajectory.Body_Set(j).ephem0.r(1)=obj.Min_Per(j)*coslat*coslong;
                    obj.Solution.Trajectory.Body_Set(j).ephem0.r(2)=obj.Min_Per(j)*coslat*sinlong;
                    obj.Solution.Trajectory.Body_Set(j).ephem0.r(3)=obj.Min_Per(j)*sinlat;
                end
            end
            
            % Now Update the Trajectory
            [obj.Solution, DeltaV ]= obj.Solution.Compute_DeltaV(tin(1:obj.Solution.Trajectory.Nbody));
            gradient=[];
            
            % Remember to update the value of tlast1!!!
            tlast1=tin;
            
            % The best permutation of Transfers is 'Best'
            Best =  obj.Solution.Trajectory.Best;
            
            % IF the value of C3 is constrained at the home planet:
            for i=1:obj.Solution.Trajectory.Nbody
                if (obj.Max_dV(i)<1e50)
                    dVeqold(i)= obj.Solution.Trajectory.dV(Best,i) - obj.Max_dV(i);
                end
            end
          
            % Now calculate minimum periapsis constraints for each body not
            % at beginning or end of Traejctory: Remember to ignore if
            % There is no Encounter for whatever reason.
            for j = 2:obj.Solution.Trajectory.Ntrans
                if(bitand(obj.Solution.Trajectory.NO_ENCOUNTER,2^j))
                    ceqold(j)=0;
                else
                    ceqold(j) = obj.Solution.Trajectory.Hyperbola(Best,j).Planet.radius +obj.Min_Per(j)- obj.Solution.Trajectory.Hyperbola(Best,j).Per;
                end
            end

            % Now calculate Maximum Periapsis Constraints for each body not
            % at beginning or end of Traejctory: Remember to ignore if
            % There is no Encounter for whatever reason.
            for j=(obj.Solution.Trajectory.Ntrans+1):(2*obj.Solution.Trajectory.Ntrans-1)
               pointer = j-obj.Solution.Trajectory.Ntrans+1;
               obj.Solution.Trajectory.Body_Set(pointer)=obj.Solution.Trajectory.Body_Set(pointer).Sphere_Of_Influence();
              SphereI = obj.Solution.Trajectory.Body_Set(pointer).SpoI;

                if(bitand(obj.Solution.Trajectory.NO_ENCOUNTER,2^pointer))
                    ceqold(j)=0;
                else
                    ceqold(j) =obj.Solution.Trajectory.Hyperbola(Best,pointer).Per- SphereI;
                end

            end
            
            % Finally calculate Minimum Perihelia Constraints where present
            for j = 1:obj.Solution.Trajectory.Ntrans
               if (bitand(obj.Perihelia_flag,2^j))
                   PER= obj.Solution.Trajectory.Trans_Set(j).transfer_body(obj.Solution.Trajectory.perm(Best,j)).orbit.a*(1-obj.Solution.Trajectory.Trans_Set(j).transfer_body(obj.Solution.Trajectory.perm(Best,j)).orbit.e);
                   reqold(j) = obj.Perihelia(j+1)-PER;
               else
                   obj.Solution.Trajectory.Trans_Set(j)=obj.Solution.Trajectory.Trans_Set(j).Calculate_Perihelion();
                   reqold(j) = obj.Perihelia(j+1)-obj.Solution.Trajectory.Trans_Set(j).perihelion(obj.Solution.Trajectory.perm(Best,j));
               end
            end
            
            
            return;
        end

        
        % Periapsis Constraints
        
        function [  con, gradient ] = Per_NLopt(tin,run_mode)
%# Per_NLopt                :   Calculates Periapsis Constraints for Optimize_Mission
        if ~isequal(tin,tlast1)
                [DeltaV,gradient] = Update_Traj(tin);
        end
           % Periapsis Constraints
           ceq=ceqold;
           con=ceq(run_mode);
           gradient = [];
       return;
        end
        
        
         % DeltaV Constraints
        
        function [  con, gradient ] = dV_NLopt(tin,run_mode)
%# dV_NLopt                :   Calculates DeltaV Constraints for Optimize_Mission
        if ~isequal(tin,tlast1)
                [DeltaV,gradient] = Update_Traj(tin);
        end
           % Periapsis Constraints
           dVeq=dVeqold;
           con=dVeq(run_mode);
           gradient = [];
       return;
        end      
        
        
        % Perihelion Constraints

       
        function [  con, gradient ] = Perhel(tin,run_mode)
 %# Perhel                   :   Calculates Perihelion Constraints for Optimize_Mission 
 
            if ~isequal(tin,tlast1)
                    [DeltaV,gradient] = Update_Traj(tin);
            end
            
            req=reqold;
                
            % Perihelion Constraints
           con=req(obj.Per_Pointer(run_mode)-1); 
          % con/obj.AU
           gradient = [];
       return;
        end

       
        function [cond,gradient] = Overall_Duration(tin)
%# Overall_Duration         :   Calculates Overall Mission Duration Constraints for Optimize_Mission 
                
                cond = 0;
                for j = 2:obj.Solution.Trajectory.Nbody
                    cond=cond+tin(j);
                end
                cond = cond-obj.Max_Duration;
             
            gradient = [];
            return;
        end
    end
    
    % Output Results
    

    function obj = View_Results(obj, numdata, Runmode)
%# View_Results             :   Displays Interplanetary Trajectory in 2D & 3D Plot form as well as orbits
%#                          :   of SSO's. Also Displays Speed Against Radial Distance Plot of Trajectory

    if (Runmode==1)
       PlotMiss = obj.Global_Solution.Trajectory;
    else
       PlotMiss = obj.Current_Mission.Trajectory;
    end
    

    % Plot Transfers
    
    nplanets =PlotMiss.Nbody;
    ntrans =PlotMiss.Ntrans;

    % Specify Time Range

    X=zeros(nplanets,numdata,3);

%     Firstly Planets
    for i=1:nplanets
        if PlotMiss.Body_Set(i).orbit.e>=1
            Time_Range =  PlotMiss.Trans_Set(i-1).tar-PlotMiss.Trans_Set(i-1).td ;
        else
            Time_Range = PlotMiss.Body_Set(i).orbit.TP;
            
        end
        if i==1
            Time_Range=min(Time_Range,-PlotMiss.Trans_Set(i).td+obj.Max_Spice_Select(i));
            tplot=linspace(PlotMiss.Trans_Set(i).td,PlotMiss.Trans_Set(i).td+Time_Range,numdata);
        else
            Time_Range=min(Time_Range,-PlotMiss.Trans_Set(i-1).td+obj.Max_Spice_Select(i));
            tplot=linspace(PlotMiss.Trans_Set(i-1).td,PlotMiss.Trans_Set(i-1).td+Time_Range,numdata);
        end
        for j=1:numdata
          if (bitand(2^i,obj.Current_Mission.Out_Of_Spice_Bounds))
              mode=1;
          else
              mode=2;
          end
         
          PlotMiss.Body_Set(i)=PlotMiss.Body_Set(i).compute_ephem_at_t(tplot(j),mode,1e-4);
            X(i,j,1)=PlotMiss.Body_Set(i).ephemt.r(1);
            X(i,j,2)=PlotMiss.Body_Set(i).ephemt.r(2);
            X(i,j,3)=PlotMiss.Body_Set(i).ephemt.r(3);
          %  if i==2
              %  tplot(j)
          %      180/pi*acos(dot(PlotMiss.Body_Set(i).ephemt.r,PlotMiss.Trans_Set(3).ephema(1).r)/PlotMiss.Body_Set(i).ephemt.R/PlotMiss.Trans_Set(3).ephema(1).R)
          %  end
        end
    end

    Y1=zeros(ntrans,numdata);
    Y2=zeros(ntrans,numdata);
    Y3=zeros(ntrans,numdata); 
    
    HT=zeros(1,numdata*ntrans);
    HR=zeros(1,numdata*ntrans);
    HV=zeros(1,numdata*ntrans);
   
    % Secondly Transfers
    
    for i=1:ntrans
        Best_Perm =PlotMiss.perm(PlotMiss.Best,i);
        tt=linspace(PlotMiss.Trans_Set(i).td,PlotMiss.Trans_Set(i).tar,numdata);
      
        for j=1:numdata
          PlotMiss.Trans_Set(i).transfer_body(Best_Perm)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).compute_ephem_at_t(tt(j),1,1e-10);

                HT((i-1)*numdata+j)=tt(j);
                HR((i-1)*numdata+j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.R/obj.AU;
                HV((i-1)*numdata+j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.V;

                if((j==1)&&(i>1))
                    HV((i-1)*numdata+j)=HV((i-1)*numdata+j)+PlotMiss.dV(PlotMiss.Best,i);
                end
            Y1(i,j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.r(1);
            Y2(i,j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.r(2);
            Y3(i,j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.r(3);
         %   Time=datetime(tt(j)/24/60/60,'ConvertFrom','juliandate');
        %    Disp = sprintf('%d Time= %s Y1=%f Y2=%f\n' ,j,Time,Y1(i,j),Y2(i,j));
            
        %    disp(Disp);
        end
    end
   %PlotMiss.Trans_Set(i).transfer_body.ephemt
   %PlotMiss.Body_Set(nplanets).ephemt

 figure(1);           
 axis equal;
    for i=1:nplanets
        plot(-X(i,:,2),X(i,:,1),'--');
        hold on;
    end
    for i=1:ntrans
        plot(-Y2(i,1:numdata),Y1(i,1:numdata));
        hold on;
    end
    hold off;
    figure(2);
    plot(HR,HV);
    
    figure(3);
    axis equal;
    for i=1:nplanets
        plot3(-X(i,:,2),X(i,:,1),X(i,:,3),'--');
        hold on;
    end
    for i=1:ntrans
        plot3(-Y2(i,1:numdata),Y1(i,1:numdata),Y3(i,1:numdata));
        hold on;
    end
    hold off;
   return;

    end


function obj = View_Info(obj,Runmode)
%# View_Info                :   Display Basic Info in numeric form about the encounter with each SSO in turn

    if (Runmode==1)
         PrinMiss = obj.Global_Solution;
    else
         PrinMiss = obj.Current_Mission;
    end
    
    f=figure(10);
    f.Position = [ 50 50 1650 700 ];
    f.Name = 'Velocity and DeltaV Information Together with Periapsis for each Solar System Object Visited';
    
    D{1} = "";
    D{2} = "";
    D{3} = "   Number      Planet                 Time        Arrival speed      Departure speed   DeltaV     Cumulative DeltaV      Periapsis";
    D{4} = "                                                      m/s                 m/s            m/s           m/s                  km   ";
    D{5} = ""; 
    
    cumdV =0;
    for i = 1:PrinMiss.Trajectory.Nbody
        Time = cspice_et2utc(PrinMiss.Absolute_Times(i),'C',0);
        if i==1||(bitand(PrinMiss.Trajectory.NO_ENCOUNTER,2^i))
            index=1;
            Periapsis = "N/A";
        elseif (bitand(PrinMiss.Trajectory.NO_ENCOUNTER,2^i))
            index=PrinMiss.Trajectory.perm(PrinMiss.Trajectory.Best,i-1);
            Periapsis = "N/A";
        else
            index=PrinMiss.Trajectory.perm(PrinMiss.Trajectory.Best,i-1);
            
        end
        if i==PrinMiss.Trajectory.Nbody
            index2=PrinMiss.Trajectory.perm(PrinMiss.Trajectory.Best,i-1);
            if (PrinMiss.FlybyRendez==1)
               DELTAV = PrinMiss.Trajectory.dV(PrinMiss.Trajectory.Best,i);
            else
                DELTAV = 0;
            end
            cumdV = cumdV + DELTAV;
            Periapsis = "N/A";
        else
            index2=PrinMiss.Trajectory.perm(PrinMiss.Trajectory.Best,i);
            DELTAV = PrinMiss.Trajectory.dV(PrinMiss.Trajectory.Best,i);
            cumdV = cumdV + DELTAV;
            if(i>1&&~bitand(PrinMiss.Trajectory.NO_ENCOUNTER,2^i))
                Periapsis = sprintf("%8.1f",(PrinMiss.Trajectory.Hyperbola(PrinMiss.Trajectory.Best,i).Per- PrinMiss.Trajectory.Body_Set(i).radius)/1000.0);
            end
         end
        
        D{i+5} = sprintf("     %d %20s %22s     %8.1f       %8.1f      %8.1f         %8.1f           %10s",i,PrinMiss.Trajectory.Body_Set(i).name,Time,norm(PrinMiss.Trajectory.VA(:,i,index)),norm(PrinMiss.Trajectory.VD(:,i,index2)),DELTAV,cumdV,Periapsis);
        

        
    %    disp(D);
    %    disp(E)
    end
    for j=PrinMiss.Trajectory.Nbody+6:obj.Max_NBody+6
        D{j}="";
    end
    u=uicontrol('Style','edit','Min',1,'Max',3);
    u.FontName = 'Courier';
    u.FontSize = 14;
    u.HorizontalAlignment = 'left';
    
   [outstring, ~]  = textwrap( u, D , 200);

   set(u,'String',outstring, 'Position', [10 10 1600 675]);

end    

   
function obj = View_Planetary_Encounters(obj, numdata,Runmode )
%# View_Planetary_Encounters:   Displays Plots of Distance Against Time for Each SSO Encounter as well
%#                          :   As 3D Trajectory for each SSO Encounter    

    if (Runmode==1)
         EncMiss = obj.Global_Solution.Trajectory;
         TCLOSEST = obj.Global_Solution.Absolute_Times;
    else
         EncMiss = obj.Current_Mission.Trajectory;
         TCLOSEST = obj.Current_Mission.Absolute_Times;
    end
       % Plot Transfers

    ntrans =EncMiss.Ntrans;

    % Specify Time Range

    XI=zeros(numdata,3);
    XO=zeros(numdata,3);
    RI=zeros(1,numdata);
    RO=zeros(1,numdata);
    RMAX=zeros(1,ntrans-1);
    TSTART=zeros(1,ntrans-1);
    TEND=zeros(1,ntrans-1);
    
    for i=2:ntrans
        if(bitand(EncMiss.NO_ENCOUNTER,2^i))
            continue;
        end
        j=i-1;
        Best_Perm =EncMiss.perm(EncMiss.Best,i);
        EncMiss.Hyperbola(Best_Perm,i).Planet = EncMiss.Hyperbola(Best_Perm,i).Planet.Sphere_Of_Influence();
        EncMiss.Hyperbola(Best_Perm,i) = EncMiss.Hyperbola(Best_Perm,i).Orbits_From_Hyperbolas();
        
        % RMAX is set to Laplace Sphere of Influence
        
        RMAX(j) = EncMiss.Hyperbola(Best_Perm,i).Planet.SpoI/100;
        
        
        % HYPS is set to Hyperbola
        
        HYPS = EncMiss.Hyperbola(Best_Perm,i);
        
        % PROBIN is set to Arriving Probe before Periapsis
        
        PROBIN = HYPS.Probe;
        PROBIN.orbit;
        
        % PROBOU is set to Departing Probe after Periapsis
        
        PROBOU = HYPS.Probe2 ;
        
        % Initialize Probes
        
        PROBIN = PROBIN.orbittoephem(0) ;       
        PROBOU = PROBOU.orbittoephem(0);

        % Start Time to Encounter Plot
        
        TSTART(j) = -RMAX(j)/norm(HYPS.VA);
        
        % Finish Time to Encounter Plot
        
        TEND(j) = -TSTART(j);
        
        tt = linspace(TSTART(j),0,numdata);
        
        
        datei = cspice_et2utc(tt+TCLOSEST(i),'C',0); 
        
        for k=1:numdata
            PROBIN=PROBIN.compute_ephem_at_t( tt(k), 1, 1e-10);
            XI(k,1)=PROBIN.ephemt.r(1,1);
            XI(k,2)=PROBIN.ephemt.r(2,1);
            XI(k,3)=PROBIN.ephemt.r(3,1);
            RI(k)=PROBIN.ephemt.R;
        end

        tt2 = linspace(0,TEND(j),numdata);
        
        dateo = cspice_et2utc(tt2+TCLOSEST(i),'C',0);
        
        for k=1:numdata
            PROBOU=PROBOU.compute_ephem_at_t( tt2(k), 1, 1e-10);
            XO(k,1)=PROBOU.ephemt.r(1,1);
            XO(k,2)=PROBOU.ephemt.r(2,1);
            XO(k,3)=PROBOU.ephemt.r(3,1);
            RO(k)=PROBOU.ephemt.R;
        end
        
        figure(30+j);
        titlestring=sprintf('Fly-by of Planet %d is %s Periapsis on %s', i, EncMiss.Body_Set(i).name,cspice_et2utc(TCLOSEST(i),'C',0));
       
        plot(datetime(datei,'InputFormat','yyyy MMM dd HH:mm:ss'),RI);
        hold on;
        plot(datetime(dateo,'InputFormat','yyyy MMM dd HH:mm:ss'),RO);
        title(titlestring);
        hold off;
        figure(10+j);
        plot3(XI(:,1),XI(:,2),XI(:,3));
        hold on;
        plot3(XO(:,1),XO(:,2),XO(:,3));
        title(titlestring);
        hold off;
    end

end


function obj = View_DeltaV_Vs_Time(obj,numdata,Runmode,TRANGE)
%# View_DeltaV_Vs_Time      :   Displays Plots of DeltaV Against Time

    if (Runmode==1)
         DVMiss = obj.Global_Solution;
    else
         DVMiss = obj.Current_Mission;
    end
    
  %  tt = linspace(DVMiss.Mission_Times(1)-TRANGE/2,DVMiss.Mission_Times(1)+TRANGE/2,numdata);
  %  datev = strings(numdata);
  %  obj.Local_Solution = DVMiss;
  %  for i=1:numdata
  %      obj.Local_Solution.Mission_Times(1) = tt(i);
  %      datev{i} = cspice_et2utc(tt(i),'C',0);
  %      obj.factor=0.5;
  %      obj = obj.Optimize_Mission(2);
  %      TotaldV(i)=obj.Local_Solution.TotaldV;
  %  end
     
  %  tt = linspace(obj.Min_time(1),obj.Max_time(1),numdata);
  %  tt2 = linspace(730*24*60*60,4000*24*60*60,numdata);
  %  datev = strings(numdata,numdata);
  %  tt3 = zeros(numdata,numdata);
  %  TotaldV= zeros(numdata,numdata);
    
    tt = linspace(obj.Min_time(1),obj.Max_time(1),numdata);
    tt2 = linspace(730*24*60*60,4000*24*60*60,numdata);
    datev = strings(numdata);
    tt3 = zeros(numdata);
    TotaldV= zeros(numdata);
        
   % for i=1:numdata
        
  %      for j=1:numdata
  %          datev{i,j} = cspice_et2utc(tt(i),'C',0);
  %          obj.Current_Mission.Mission_Times(1) = tt(i);
  %          obj.Current_Mission.Mission_Times(2) = tt2(j);
  %          tt3(i,j)=tt2(j);
  %          [obj.Current_Mission TotaldV(i,j)] = obj.Current_Mission.Compute_DeltaV( obj.Current_Mission.Mission_Times );
  %      end
  %  end 
     
  for i=1:numdata
      datev{i}=cspice_et2utc(tt(i),'C',0);
      obj.Current_Mission.Mission_Times(1) = tt(i);
      obj.Current_Mission.Mission_Times(2) = 10*24*60*60;
      obj.Min_time(1)=tt(i);
      obj.Max_time(1)=tt(i)+365*24*60*60;
      obj = obj.Optimize_Mission(4);
      [obj.Current_Mission TotaldV(i)] = obj.Current_Mission.Compute_DeltaV( obj.Current_Mission.Mission_Times );
      obj.Run_Time = 5*60;
  end
   figure(100);
  plot(datetime(datev,'InputFormat','yyyy MMM dd HH:mm:ss'),TotaldV);
 


    
  %  surf(datetime(datev,'InputFormat','yyyy MMM dd HH:mm:ss'),tt3/24/60/60,TotaldV);

 %   plot(tt3(1,:)/24/60/60,TotaldV(1,:));
 %   hold on;
 %   plot(tt3(numdata,:)/24/60/60,TotaldV(numdata,:));
 %   hold off;

end


function obj = View_Orbit_Info(obj,Runmode)
%# View_Orbit_Info          :   Displays Orbital Information in numeric form

    if (Runmode==1)
         OrbitM = obj.Global_Solution;
    else
         OrbitM = obj.Current_Mission;
    end
    f=figure(200);
    f.Position = [ 50 50 1500 700 ];
    f.Name = 'Orbital Information for Encounter of Each Solar System Object';
    
    D{1} = "";
    D{2} = "";
    D{3} = "   Number      Planet       Periapsis Time              Periapsis      Arr Ecc    Dep Ecc     Inclination         LOAN        AOP";
    D{4} = "                                                           km                                   degs             degs        degs";
    D{5} = ""; 
    
    for i = 1:OrbitM.Trajectory.Nbody
        Time = cspice_et2utc(OrbitM.Absolute_Times(i),'C',0);
        if i==1
            Data = "N/A";
        elseif(bitand(OrbitM.Trajectory.NO_ENCOUNTER,2^i))
            Data = "N/A";
        end
        if i==OrbitM.Trajectory.Nbody
            Data = "N/A";
        else
            index2=OrbitM.Trajectory.perm(OrbitM.Trajectory.Best,i);
            if(i>1&&~bitand(OrbitM.Trajectory.NO_ENCOUNTER,2^i))
                
                Periapsis = (OrbitM.Trajectory.Hyperbola(OrbitM.Trajectory.Best,i).Per- OrbitM.Trajectory.Body_Set(i).radius)/1000;
                OrbitM.Trajectory.Hyperbola(index2,i) = OrbitM.Trajectory.Hyperbola(index2,i).Orbits_From_Hyperbolas();
                EccA = OrbitM.Trajectory.Hyperbola(index2,i).Probe.orbit.e;
                EccD = OrbitM.Trajectory.Hyperbola(index2,i).Probe2.orbit.e;
                Inc = OrbitM.Trajectory.Hyperbola(index2,i).Probe.orbit.I * 180/pi;
                LOAN = OrbitM.Trajectory.Hyperbola(index2,i).Probe.orbit.loan * 180/pi;
                AOP = OrbitM.Trajectory.Hyperbola(index2,i).Probe.orbit.aop * 180/pi;
                Data = sprintf("%12.1f     %6.3f     %6.3f        %7.2f        %7.2f    %7.2f", Periapsis,EccA,EccD,Inc,LOAN,AOP);
            end
         end
        
        D{i+5} = sprintf("     %d %20s %22s     %s",i,OrbitM.Trajectory.Body_Set(i).name,Time,Data);
  
    end
    for j=OrbitM.Trajectory.Nbody+6:obj.Max_NBody+6
        D{j}="";
    end
    u=uicontrol('Style','edit','Min',1,'Max',3);
    u.FontName = 'Courier';
    u.FontSize = 14;
    u.HorizontalAlignment = 'left';
    
   [outstring, ~]  = textwrap( u, D , 200);

   set(u,'String',outstring, 'Position', [10 10 1600 675]);
    
        f=figure(300);
    f.Position = [ 50 50 1650 700 ];
    f.Name = 'Interplanetary Transfer Orbital Information';
    
    E{1} = "";
    E{2} = "";
    E{3} = "   Transfer    Dep Time             Arr Time             PerihelionR   PerihelionO    Ecc        Inclination      LOAN        AOP";
    E{4} = "                                                              AU            AU                        degs        degs        degs";
    E{5} = ""; 
    
    for i = 1:OrbitM.Trajectory.Ntrans
        Time1 = cspice_et2utc(OrbitM.Absolute_Times(i),'C',0);
        Time2 = cspice_et2utc(OrbitM.Absolute_Times(i+1),'C',0);

            index2=OrbitM.Trajectory.perm(OrbitM.Trajectory.Best,i);
                OrbitM.Trajectory.Trans_Set(i)=OrbitM.Trajectory.Trans_Set(i).Calculate_Perihelion();
                PerihelionR = OrbitM.Trajectory.Trans_Set(i).perihelion(index2)/obj.AU;
                 
                Ecc = OrbitM.Trajectory.Trans_Set(i).transfer_body(index2).orbit.e;
                PerihelionO =OrbitM.Trajectory.Trans_Set(i).transfer_body(index2).orbit.a*(1-Ecc)/obj.AU;
                Inc = OrbitM.Trajectory.Trans_Set(i).transfer_body(index2).orbit.I*180/pi;
                LOAN = OrbitM.Trajectory.Trans_Set(i).transfer_body(index2).orbit.loan*180/pi;
                AOP = OrbitM.Trajectory.Trans_Set(i).transfer_body(index2).orbit.aop*180/pi;
                Data = sprintf("%10.4f  %10.4f    %6.3f         %5.2f          %5.2f      %5.2f", PerihelionR,PerihelionO,Ecc,Inc,LOAN,AOP);
        
        E{i+5} = sprintf("  %d - %d | %20s | %20s |   %s",i,i+1,Time1,Time2,Data);
  
    end
    
    for j=OrbitM.Trajectory.Nbody+6:obj.Max_NBody+6
        E{j}="";
    end
    
    u=uicontrol('Style','edit','Min',1,'Max',3);
    u.FontName = 'Courier';
    u.FontSize = 14;
    u.HorizontalAlignment = 'left';
    
   [outstring, ~]  = textwrap( u, E , 150);
    
   set(u,'String',outstring, 'Position', [10 10 1600 650]);
   
end    


    function obj = Animate_Results(obj, numdata, Runmode,titleanim)
%# Animate_Results          :   Animates Interplanetary Trajectory

    if (Runmode==1)
       PlotMiss = obj.Global_Solution.Trajectory;
    else
       PlotMiss = obj.Current_Mission.Trajectory;
    end
    

    % Plot Transfers
    
    nplanets =PlotMiss.Nbody;
    ntrans =PlotMiss.Ntrans;

    % Specify Time Range

    X=zeros(nplanets,numdata,3);

%     Firstly Planets
    for i=1:nplanets
        
       if (PlotMiss.Body_Set(i).orbit.e>=1)||isempty(PlotMiss.Body_Set(i).orbit.e)
           TIME=  PlotMiss.Trans_Set(i-1).tar-PlotMiss.Trans_Set(i-1).td; 
       else
           TIME=PlotMiss.Body_Set(i).orbit.TP;
       end

       if i==1
            tplot=linspace(PlotMiss.Trans_Set(i).td,PlotMiss.Trans_Set(i).td+TIME,numdata);
       else
            tplot=linspace(PlotMiss.Trans_Set(i-1).td,PlotMiss.Trans_Set(i-1).td+TIME,numdata);
       end
       
       for j=1:numdata
            if (bitand(2^i,obj.Current_Mission.Out_Of_Spice_Bounds))
                mode=1;
            else
                mode=2;
            end
          PlotMiss.Body_Set(i)=PlotMiss.Body_Set(i).compute_ephem_at_t(tplot(j),mode,1e-4);
            X(i,j,1)=PlotMiss.Body_Set(i).ephemt.r(1);
            X(i,j,2)=PlotMiss.Body_Set(i).ephemt.r(2);
            X(i,j,3)=PlotMiss.Body_Set(i).ephemt.r(3);
        end
    end

    TT=zeros(ntrans,numdata);
    Y1=zeros(ntrans,numdata);
    Y2=zeros(ntrans,numdata);
    Y3=zeros(ntrans,numdata); 
    Xp1=zeros(ntrans,nplanets,numdata);
    Xp2=zeros(ntrans,nplanets,numdata);
    Xp3=zeros(ntrans,nplanets,numdata);
    
    HT=zeros(1,numdata*ntrans);
    HR=zeros(1,numdata*ntrans);
    HV=zeros(1,numdata*ntrans);
    HCD=zeros(1,numdata*ntrans);
   
    % Secondly Transfers
    
    CUM_DIST=0;
    for i=1:ntrans
        Best_Perm =PlotMiss.perm(PlotMiss.Best,i);
        tt=linspace(PlotMiss.Trans_Set(i).td,PlotMiss.Trans_Set(i).tar,numdata);
        
        TOLD=tt(1);
        
        for j=1:numdata
          PlotMiss.Trans_Set(i).transfer_body(Best_Perm)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).compute_ephem_at_t(tt(j),1,1e-10);
            TT(i,j)=tt(j);
            Y1(i,j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.r(1);
            Y2(i,j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.r(2);
            Y3(i,j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.r(3);
            HR(i,j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.R/obj.AU;
            HV(i,j)=PlotMiss.Trans_Set(i).transfer_body(Best_Perm).ephemt.V/1000;
            CUM_DIST=CUM_DIST+ HV(i,j)*(TT(i,j)-TOLD);
            HCD(i,j)=CUM_DIST;
            TOLD=TT(i,j);
            % Do Exact Position of planets
            for n=1:nplanets
                if (bitand(2^n,obj.Current_Mission.Out_Of_Spice_Bounds))
                    mode=1;
                else
                    mode=2;
                end
                PlotMiss.Body_Set(n)=PlotMiss.Body_Set(n).compute_ephem_at_t(tt(j),mode,1e-4);
                Xp1(i,n,j)=PlotMiss.Body_Set(n).ephemt.r(1);
                Xp2(i,n,j)=PlotMiss.Body_Set(n).ephemt.r(2);
                Xp3(i,n,j)=PlotMiss.Body_Set(n).ephemt.r(3);
            end
                
        end
    end


   MAXFLIGHTTIME= max(obj.Current_Mission.Mission_Times(2:PlotMiss.Nbody));
   
   MAXY1=max(max(Y1));
   MINY1=min(min(Y1));
   
   MAXY2=max(max(Y2));
   MINY2=min(min(Y2));
   
   MAXY3=max(max(Y3));
   MINY3=min(min(Y3));
   
   MAXY1ABS= max(abs(MAXY1),abs(MINY1));
   MAXY2ABS= max(abs(MAXY2),abs(MINY2));
   MAXY3ABS= max(abs(MAXY3),abs(MINY3));
   MAXTOT=max(MAXY1ABS,MAXY2ABS);
   
 
    figure('Position', [0 0 900 900]);
    titlestr='';
    for i=1:nplanets
        titlestr=sprintf('%s %s',titlestr,PlotMiss.Body_Set(i).name);
    end
    

    for i=1:nplanets
        axis([-1.5*MAXTOT 1.5*MAXTOT -1.5*MAXTOT 1.5*MAXTOT]);
  %      axis([-0.1*MAXTOT 0.1*MAXTOT -0.1*MAXTOT 0.1*MAXTOT]);
        ax=gca;
        ax.Position=[0 0 1 1];
%        axis([-0.25*MAXTOT 0.25*MAXTOT -0.25*MAXTOT 0.25*MAXTOT -1e12 1e12]);
%        axis([-1.5*MAXTOT 1.5*MAXTOT -1.5*MAXTOT 1.5*MAXTOT -1.5*MAXTOT 1.5*MAXTOT]);
        plot(-X(i,:,2),X(i,:,1),':');
%         plot3(-X(i,:,2),X(i,:,1),X(i,:,3),'--');
        hold on;
    end
    TIT=title({' ';titleanim});
    set(TIT,'VerticalAlignment','top','HorizontalAlignment', 'center');
    drawnow;
    
    myVideo=VideoWriter('TrajVideo.mp4','MPEG-4');
    open(myVideo);
    kcum=0;
    CUMDV=0;
    
    k=zeros(1,ntrans);
    p=plot(1:nplanets+1);

    for i=1:ntrans
        interval= int64(MAXFLIGHTTIME/obj.Current_Mission.Mission_Times(i+1));
        
        datestr=cspice_et2utc(TT(i,1),'C',0);
        planetstr = PlotMiss.Body_Set(i).name;
        legendstr= sprintf('%s %s',planetstr,datestr);
        
        CUMDV=CUMDV+PlotMiss.dV(PlotMiss.Best,i);
        
        for j=1:interval:numdata
            k(i)=k(i)+1;
            
            if (i>1)
                update = sprintf('%s\nCLOSEST APPROACH = %11.2fkm',legendstr,(PlotMiss.Hyperbola(PlotMiss.Best,i).Per-PlotMiss.Body_Set(i).radius)/1000);
                an1=annotation('textbox',[.4 .12 .47 .04],'String',update,'FontName','FixedWidth','LineStyle','none');
            else
                an1=annotation('textbox',[.4 .12 .47 .04],'String',legendstr,'FontName','FixedWidth','LineStyle','none');
            end
            CurTimestr=cspice_et2utc(TT(i,j),'C',0);
            Info1str=sprintf('DISTANCE FROM SUN = %4.1fAU\nSPEED = %4.1fkm/s',HR(i,j),HV(i,j));
            Info2str=sprintf('DISTANCE TRAVELLED = %14.0fkm',HCD(i,j));
            Info3str=sprintf('DeltaV at %s = %4.1fkm/s\nCUMULATIVE DeltaV = %4.1fkm/s',planetstr,PlotMiss.dV(PlotMiss.Best,i)/1000,CUMDV/1000);
            an2=annotation('textbox',[.4 .89 .47 .02],'String',CurTimestr,'FontName','FixedWidth','LineStyle','none');
            an3=annotation('textbox',[.4 .85 .47 .04],'String',Info1str,'FontName','FixedWidth','LineStyle','none');
            an4=annotation('textbox',[.4 .83 .47 .02],'String',Info2str,'FontName','FixedWidth','LineStyle','none');
            an5=annotation('textbox',[.4 .10 .47 .02],'String',Info3str,'FontName','FixedWidth','LineStyle','none');
            an1.Color='black';
            an2.Color='black';
            an3.Color='black';
            an4.Color='black';
            an5.Color='black';
            % Do all planets
     
            for l=1:nplanets
                if (PlotMiss.Body_Set(l).Fixed_Point>0)
                    continue;
                end
                axis([-1.5*MAXTOT 1.5*MAXTOT -1.5*MAXTOT 1.5*MAXTOT ]);
                p(l)=plot(-Xp2(i,l,j),Xp1(i,l,j),'or');
                XAN = [(0.5-Xp2(i,l,j)/3.0/MAXTOT-0.0001) (0.5-Xp2(i,l,j)/3.0/MAXTOT)];
                YAN = [(0.5+Xp1(i,l,j)/3.0/MAXTOT-0.0001) (0.5+Xp1(i,l,j)/3.0/MAXTOT)];

                a(l)=annotation('textarrow',XAN,YAN,'HeadStyle','none','String',PlotMiss.Body_Set(l).name,'Color','red','FontSize',9);
                hold on;
            end
            
            
            
            plot(-Y2(i,1:j),Y1(i,1:j),'Color','black');
            p(nplanets+1)=plot(-Y2(i,j),Y1(i,j),'o','Color','black');
        %    plot3(-Y2(i,1:j),Y1(i,1:j),Y3(i,1:j));
     %       axis([-0.1*MAXTOT 0.1*MAXTOT -0.1*MAXTOT 0.1*MAXTOT]);
       %     axis([-0.25*MAXTOT 0.25*MAXTOT -0.25*MAXTOT 0.25*MAXTOT -1e12 1e12]);
            axis([-1.5*MAXTOT 1.5*MAXTOT -1.5*MAXTOT 1.5*MAXTOT ]);
            kcum=kcum+k(i); 
            F(kcum)=getframe(gcf);
            writeVideo(myVideo, F(kcum));
            hold on;
            for l=1:nplanets
                if (PlotMiss.Body_Set(l).Fixed_Point>0)
                    continue;
                end
                %bplot(-Xp2(i,l,j),Xp1(i,l,j),'ow');
                p(l).Color='none';
                a(l).Color='none';
                hold on;
            end
            p(nplanets+1).Color='none';
            
            an1.Color='none';
            an2.Color='none';
            an3.Color='none';
            an4.Color='none';
            an5.Color='none';
        end
        hold on;
      
 
    end
    
    datestr=cspice_et2utc(TT(i,j),'C',0);
    planetstr = PlotMiss.Body_Set(i+1).name;
    legendstr= sprintf('%s %s',planetstr,datestr);
    CUMDV=CUMDV+PlotMiss.dV(PlotMiss.Best,i+1);
    CurTimestr=cspice_et2utc(TT(i,j),'C',0);
    Info1str=sprintf('DISTANCE FROM SUN = %4.1fAU\nSPEED = %4.1fkm/s',HR(i,j),HV(i,j));
    Info2str=sprintf('DISTANCE TRAVELLED = %14.0fkm',HCD(i,j));
    Info3str=sprintf('DeltaV at %s = %4.1fkm/s\nCUMULATIVE DeltaV = %4.1fkm/s',planetstr,PlotMiss.dV(PlotMiss.Best,i+1)/1000,CUMDV/1000);
    
    for k=1:80
        an1=annotation('textbox',[.4 .12 .47 .04],'String',legendstr,'FontName','FixedWidth','LineStyle','none','Color','black');
        an2=annotation('textbox',[.4 .89 .47 .02],'String',CurTimestr,'FontName','FixedWidth','LineStyle','none','Color','black');
        an3=annotation('textbox',[.4 .85 .47 .04],'String',Info1str,'FontName','FixedWidth','LineStyle','none','Color','black');
        an4=annotation('textbox',[.4 .83 .47 .02],'String',Info2str,'FontName','FixedWidth','LineStyle','none','Color','black');
        an5=annotation('textbox',[.4 .10 .47 .02],'String',Info3str,'FontName','FixedWidth','LineStyle','none','Color','black');
        kcum=kcum+1;
        F(kcum)=getframe(gcf);
        writeVideo(myVideo, F(kcum));
    end
    hold off;

    
    close(myVideo);
    'Done'
   return;

end
end
end
    
