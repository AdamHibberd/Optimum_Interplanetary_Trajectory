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
%# Per_NLopt                :   Calculates Periapsis Constraints for Optimize_Mission
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
         [row,col] = size(cover);

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
    function obj = Optimize_Mission(obj, runmode)
%# Optimize_Mission         :   Optimizes Mission provided by Current_Mission -> Solution goes to Solution           
        % Set up Inputs To Optimizer
      
        % OPTIMIZER TOOL
        if runmode == 1
            obj.Solution = obj.Current_Mission;
            % Global Minimum Selected
            opt.algorithm = NLOPT_GN_ISRES;
          %  opt.algorithm = NLOPT_GN_CRS2_LM
          %  opt.algorithm = NLOPT_GN_DIRECT_L
          %  opt.algorithm = NLOPT_GN_DIRECT_L_NOSCAL
          %  opt.algorithm = NLOPT_GN_DIRECT_L_RAND
          %  opt.algorithm = NLOPT_GN_MLSL_LDS
          %  opt.algorithm = NLOPT_GN_ORIG_DIRECT
          %  opt.algorithm = NLOPT_GN_ORIG_DIRECT_L

            opt.maxtime=obj.Run_Time;
            opt.verbose = 1;
            for i=1:obj.Solution.Trajectory.Nbody
                lb(i)=obj.Min_time(i);
                ub(i)=obj.Max_time(i);
            end
        elseif runmode == 2
            obj.Solution=obj.Local_Solution;
            % Local Minimum Selected
            opt.algorithm = NLOPT_LN_COBYLA;
        %    opt.algorithm = NLOPT_LN_AUGLAG
        %    opt.algorithm = NLOPT_LN_AUGLAG_EQ
        %    opt.algorithm = NLOPT_LN_BOBYQA
        %    opt.algorithm = NLOPT_LN_NELDERMEAD
        %    opt.algorithm = NLOPT_LN_NEWUOA_BOUND
        %    opt.algorithm = NLOPT_LN_NEWUOA
        %    opt.algorithm = NLOPT_LN_PRAXIS
        %    opt.algorithm = NLOPT_LN_SBPLX
            opt.verbose = 1;
            opt.ftol_rel = 1e-6;
        % Keep Launch Time Static in this Mode
        
            lb(1) = obj.Solution.Mission_Times(1); 
            ub(1) = obj.Solution.Mission_Times(1);
                lb(1)=obj.Min_time(1);
                ub(1)=obj.Max_time(1);  
            % Alter remaining times by factor
            obj.factor = obj.factor/2; 
            fac = obj.factor;

            for i=2:obj.Solution.Trajectory.Nbody
               lb(i)=obj.Solution.Mission_Times(i)*(1-fac);
               ub(i)=obj.Solution.Mission_Times(i)*(1+fac);
         %      obj.Min_time(i)=lb(i);
         %      obj.Max_time(i)=ub(i);
            end    
        elseif runmode == 3
            obj.Solution=obj.Current_Mission;
          opt.algorithm = NLOPT_LN_AUGLAG;
      %    opt.local_optimizer.algorithm = NLOPT_LN_COBYLA;
          opt.local_optimizer.algorithm = NLOPT_GN_CRS2_LM;
      %    opt.local_optimizer.algorithm = NLOPT_GN_DIRECT_L;
      %    opt.local_optimizer.algorithm = NLOPT_GN_DIRECT_L_NOSCAL;
       %   opt.local_optimizer.algorithm =   NLOPT_GN_DIRECT_L_RAND;
      %    opt.local_optimizer.algorithm =  NLOPT_GN_DIRECT_L_RAND_NOSCAL;
     %    opt.local_optimizer.algorithm = NLOPT_GN_DIRECT;
      %    opt.local_optimizer.algorithm =  NLOPT_GN_DIRECT_NOSCAL;
      %    opt.local_optimizer.algorithm = NLOPT_GN_ISRES;
       %   opt.local_optimizer.algorithm = NLOPT_GN_MLSL_LDS;
      %    opt.local_optimizer.algorithm = NLOPT_GN_MLSL;
      %    opt.local_optimizer.algorithm =   NLOPT_GN_ORIG_DIRECT_L;
      %    opt.local_optimizer.algorithm = NLOPT_GN_ORIG_DIRECT;
       %   opt.local_optimizer.algorithm = NLOPT_LN_BOBYQA;
      %    opt.local_optimizer.algorithm = NLOPT_LN_NELDERMEAD;
       %    opt.local_optimizer.algorithm = NLOPT_LN_NEWUOA;
         %  opt.local_optimizer.algorithm = NLOPT_LN_PRAXIS;
         % opt.local_optimizer.algorithm = NLOPT_LN_SBPLX;
         
          opt.local_optimizer.ftol_rel = 1e-6;
            opt.maxtime=obj.Run_Time;
            opt.verbose = 1;
            for i=1:obj.Solution.Trajectory.Nbody
                lb(i)=obj.Min_time(i);
                ub(i)=obj.Max_time(i);
            end
    
        
        % Keep Launch Time Static in this Mode
        
        %    lb(1) = obj.Solution.Mission_Times(1); 
        %    ub(1) = obj.Solution.Mission_Times(1);
         %       lb(1)=obj.Min_time(1);
         %       ub(1)=obj.Max_time(1);  
            % Alter remaining times by factor
         %   obj.factor = obj.factor/2;
         %   fac= obj.factor;
         %   for i=2:obj.Solution.Trajectory.Nbody
         %      lb(i)=obj.Solution.Mission_Times(i)*(1-fac);
         %      ub(i)=obj.Solution.Mission_Times(i)*(1+fac);
         %      obj.Min_time(i)=lb(i);
         %      obj.Max_time(i)=ub(i);
         %   end
        else
            obj.Solution=obj.Current_Mission;
      
            for i=1:obj.Solution.Trajectory.Nbody
                lb(i)=obj.Min_time(i);
                ub(i)=obj.Max_time(i);
            end

         %   nlcon = @(x)[ Per_NLopt(x,1); Per_NLopt(x,2) ; Per_NLopt(x,3); Per_NLopt(x,4); Per_NLopt(x,5) ];    
            
            % Set up Inequality Constraint Functions
            funcstring= '@(x)[';
            
            
            % Firstly The Minimum allowable Periapsis
            if (obj.Nconstraints>0)
                for i=1:obj.Nconstraints
                    bb_output_type{i}='PB';
                    funcstring= strcat(funcstring, sprintf(' Per_NLopt(x,%d)',i+1));
                         funcstring = strcat( funcstring ,' ;' );
                    nle(i)=-1;
                    nlrhs(i)=0;
                end
            end
 
            % Secondly The Maximum allowable Periapsis
            if (obj.Nconstraints>0)
                for i=(obj.Nconstraints+1):2*obj.Nconstraints
                    bb_output_type{i}='PB';
                    funcstring= strcat(funcstring, sprintf(' Per_NLopt(x,%d)',i+1));
                    if i<2*obj.Nconstraints
                         funcstring = strcat( funcstring ,' ;' );
                    end
                    nle(i)=-1;
                    nlrhs(i)=0;
                end
            end
            
            % Thirdly The Minimum allowable Perihelia
            if(obj.NPerihelia>0)
                funcstring = strcat( funcstring ,' ;' );
                for i=1:obj.NPerihelia
                    funcstring = strcat(funcstring,sprintf(' Perhel(x,%d)',i));
                    if i<obj.NPerihelia
                        funcstring = strcat( funcstring ,' ;' );
                    end
                    nle(i+2*obj.Nconstraints)=-1;
                    nlrhs(i+2*obj.Nconstraints)=0;
                    bb_output_type{i+2*obj.Nconstraints}='EB';
 %                   bb_output_type{i+2*obj.Nconstraints}='PB';
                end
            end
            
            Max_Dur_Flag = 0;
            % Finally the Maximum Duration
            if (obj.Max_Duration < obj.MAX_DURATION )
                Max_Dur_Flag = 1;
                funcstring = strcat( funcstring ,' ; Overall_Duration(x)');
                nle(2*obj.Nconstraints+obj.NPerihelia+1)=-1;
                nlrhs(2*obj.Nconstraints+obj.NPerihelia+1)=0;
                bb_output_type{2*obj.Nconstraints+obj.NPerihelia+1}='EB';
            end
            
            % Check to see if Limit on C3 At Home Planet
            if (obj.Min_Per(1) > 0.0 && obj.Solution.Trajectory.Body_Set(1).Fixed_Point == 0)
                funcstring = strcat(funcstring,sprintf(' ; Per_NLopt(x,%d)',1));
                nle(2*obj.Nconstraints+obj.NPerihelia+Max_Dur_Flag+1)=-1;
                nlrhs(2*obj.Nconstraints+obj.NPerihelia+Max_Dur_Flag+1)=0;
                bb_output_type{2*obj.Nconstraints+obj.NPerihelia+Max_Dur_Flag+1}='PB';
            end
            % Construct Function Handle From String.
           funcstring=strcat(funcstring, ' ]');
           nlcon = eval(funcstring)

          %  nlcon = {@(x) (Per_NLopt(x,1)) ;@(x) (Per_NLopt(x,2));@(x) (Per_NLopt(x,3));@(x) (Per_NLopt(x,4));@(x) (Per_NLopt(x,5))}
            optiSolver('NLP');
          %  nopts = nomadset('direction_type','ORTHO 2N','vns_search',0.75,'bb_output_type',bb_output_type );
         %  nopts = nomadset('direction_type','ORTHO 2N','vns_search',0.75 );
         
            % nopts = nomadset('bb_output_type',bb_output_type );
            % opts=optiset('solver','nomad','display','iter','solverOpts',nopts);          
            
            % Initial Guess           
  
             % Check for Presence of Intermediate Point
             NIP=-1;
             for i=1:obj.Solution.Trajectory.Nbody
                 if obj.Solution.Trajectory.Body_Set(i).Fixed_Point==1
                     NIP=NIP+2;
                     if NIP==1
                        tin = cat ( 2, obj.Solution.Mission_Times ,[ 0 0 ] );
                     else
                        tin = cat ( 2, tin , [ pi*i/20 0 ] );
                     end
                     lb(obj.Solution.Trajectory.Nbody+NIP)=0;
                     lb(obj.Solution.Trajectory.Nbody+NIP+1)=-pi/2;
                     ub(obj.Solution.Trajectory.Nbody+NIP)=2*pi;
                     ub(obj.Solution.Trajectory.Nbody+NIP+1)=pi/2;                     
                 end
             end
 
            for i=1:obj.Solution.Trajectory.Nbody
                tin(i) = obj.Solution.Mission_Times(i);
            end
            
            tlast1 = tin - tin;
            tlast2 = tlast1;
            tlast3 = tlast2;
            DeltaVold=0;
            condold = 0;
            ceq=zeros(1,min(1,2*obj.Nconstraints));
            ceqold=zeros(1,min(1,2*obj.Nconstraints));
            req=zeros(1,obj.Solution.Trajectory.Nbody-1);
            reqold=zeros(1,obj.Solution.Trajectory.Nbody-1); 
            VIOLATION=0;
        
            % Initialise Trajecory
            DeltaV =  Compute_DeltaV_NLopt(tin)            
          %  obj.Current_Mission = obj.Solution;
          %  obj.Global_Solution = obj.Solution;
          %   return;
            % Initialise NOMAD Optimising Settings
            if (obj.Nconstraints>0||obj.NPerihelia>0)
         %       nopts = nomadset('bb_output_type',bb_output_type ,'vns_search',0.95,'max_eval',60000)
                nopts = nomadset('bb_output_type',bb_output_type ,'vns_search',0.75,'max_eval',60000)
            elseif obj.Max_Duration < obj.MAX_DURATION
                nopts = nomadset('bb_output_type',bb_output_type,'vns_search',0.75);
            else
                nopts=nomadset('vns_search',0.75);
                nlrhs=[];
                nle=[];
            end
      %      return;
      %      opts=optiset('solver','nomad','display','iter','maxfeval',60000,'maxtime',obj.Run_Time,'solverOpts',nopts,'tolrfun',0.005);
            opts=optiset('solver','nomad','display','iter','maxfeval',60000,'maxtime',obj.Run_Time,'solverOpts',nopts); 
            Opt = opti('fun',@Compute_DeltaV_NLopt,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'x0',tin,'options',opts)
            % Run Optimization
            [Optimum,DeltaV,exitflag,info] = solve(Opt,tin) 
          % tin
          %  Optimum
            tlast1 = tin - tin;
            
            DeltaV          
            
            
            %Set Up the Trajectory for storing
            for i=1:obj.Solution.Trajectory.Nbody
                obj.Solution.Mission_Times(i) = Optimum(i);
            end
            DeltaV =  Compute_DeltaV_NLopt(Optimum)
        
            % Set Solution 
            obj.Current_Mission = obj.Solution;
            obj.Global_Solution = obj.Solution;
return;
        end
        
        % First of all Lower and Upper Bounds for Times
        opt.lower_bounds=lb;
        opt.upper_bounds=ub;
        
        % Set up Objective function i.e. minimum Total DeltaV
        opt.min_objective = @Compute_DeltaV_NLopt;
        
        % Set up inequality constraints and tolerances
        if obj.Nconstraints>0
            for i=1:obj.Nconstraints
                opt.fc{i} = (@(x) Per_NLopt(x,i));
                obj.Constr_Tol(i)= 1e-10;
                opt.fc_tol(i) = obj.Constr_Tol(i);
            end
        end        

        
        % Set up Constraint of Overall Duration of Mission
        if (obj.Max_Duration < obj.MAX_DURATION )
            opt.fc{obj.Nconstraints+obj.NPerihelia+1} = (@(x) Overall_Duration(x));
            obj.Constr_Tol(obj.Nconstraints+obj.NPerihelia+1) = 1e-10 ;
            opt.fc_tol(obj.Nconstraints+obj.NPerihelia+1)= obj.Constr_Tol(obj.Nconstraints+obj.NPerihelia+1);
        end
        
        
        
        
        % Initial Guess
        tin = obj.Solution.Mission_Times;
        tlast1 = tin - tin;
        tlast2 = tlast1;
        tlast3 = tlast2;
        DeltaVold=0;
        condold = 0;
        ceq=zeros(1,min(1,obj.Nconstraints));
        ceqold=zeros(1,min(1,obj.Nconstraints));
        req=zeros(1,obj.Solution.Trajectory.Nbody-1);
        reqold=zeros(1,obj.Solution.Trajectory.Nbody-1);       
        %Initialise Trajecory

      
         DeltaV =  Compute_DeltaV_NLopt(tin);
         
         % Coverge to Solution




    

        [Optimum, DeltaV, retcode] = nlopt_optimize(opt,tin);
      %  DeltaV
        %Set Up the Trajectory for storing
        obj.Solution.Mission_Times = Optimum;
        DeltaV =  Compute_DeltaV_NLopt(Optimum)
        if obj.Nconstraints > 0
            for i = 1:obj.Nconstraints
                [PERIAPSIS GRAD] = Per_NLopt(Optimum,i);
                PERIAPSIS
            end
        end
        % Set Solution 
        
        if (runmode ==1 || runmode ==3 )
            obj.Current_Mission = obj.Solution;
        else
            obj.Local_Solution = obj.Solution;
        end
        % Objective Function is Total DeltaV
        
        function [DeltaV, gradient  ] = Compute_DeltaV_NLopt(tin1)
%# Compute_DeltaV_NLopt     :   Calculates DeltaV as determined by times input by Optimize_Mission       
                if ~isequal(tin1,tlast1)
                    NIP=-1;
                    for j=1:obj.Solution.Trajectory.Nbody
                        if obj.Solution.Trajectory.Body_Set(j).Fixed_Point==1
                            NIP=NIP+2;
                            coslong=cos(tin1(obj.Solution.Trajectory.Nbody+NIP));
                            sinlong=sin(tin1(obj.Solution.Trajectory.Nbody+NIP));
                            coslat=cos(tin1(obj.Solution.Trajectory.Nbody+NIP+1));
                            sinlat=sin(tin1(obj.Solution.Trajectory.Nbody+NIP+1));
                            obj.Solution.Trajectory.Body_Set(j).ephem0.r(1)=obj.Min_Per(j)*coslat*coslong;
                            obj.Solution.Trajectory.Body_Set(j).ephem0.r(2)=obj.Min_Per(j)*coslat*sinlong;
                            obj.Solution.Trajectory.Body_Set(j).ephem0.r(3)=obj.Min_Per(j)*sinlat;
                        end
                    end
                    [obj.Solution, DeltaV ]= obj.Solution.Compute_DeltaV(tin1(1:obj.Solution.Trajectory.Nbody));
                   % if obj.Solution.Trajectory.DIVERGING>0
                   %     'DIVERGING'
                   % end
                   % if obj.Solution.Trajectory.NO_ENCOUNTER>0
                   %      'NO ENCOUNTER DYNAMICS for Body:'
                   %         obj.Solution.Trajectory.NO_ENCOUNTER 
                   % end
                    DeltaVold=DeltaV;
               
                    Best =  obj.Solution.Trajectory.Best;
                    VIOLATION=0;
                    if (obj.Min_Per(1)>0.0)
                        ceq(1)= obj.Solution.Trajectory.dV(Best,1)/1000 - sqrt(obj.Min_Per(1));
                    end
                    for j = 2:obj.Solution.Trajectory.Ntrans
                        if(bitand(obj.Solution.Trajectory.NO_ENCOUNTER,2^j))
                            ceq(j)=0;
                        else
                            ceq(j) = obj.Solution.Trajectory.Hyperbola(Best,j).Planet.radius +obj.Min_Per(j)- obj.Solution.Trajectory.Hyperbola(Best,j).Per;

           %                 if ceq(j) > 0 
            %                    DeltaV = 1e50*ceq(j)/(obj.Solution.Trajectory.Hyperbola(Best,j).Planet.radius +obj.Min_Per(j));
             %                   DeltaVold=DeltaV;
              %              end
                        end
             %           obj.Solution.Trajectory.Body_Set(j)=obj.Solution.Trajectory.Body_Set(j).Sphere_Of_Influence();
             %           SphereI = obj.Solution.Trajectory.Body_Set(j).SpoI;
             %           square1=((obj.Solution.Trajectory.Hyperbola(Best,j).Planet.radius +obj.Min_Per(j))/SphereI)^2;
             %           square2=(obj.Solution.Trajectory.Hyperbola(Best,j).Per/SphereI)^2;
             %           ceq(j)=square1-sign(obj.Solution.Trajectory.Hyperbola(Best,j).Per)*square2;
             %           if ceq(j) >0
             %              VIOLATION=1;
             %           end
                    end
                   
                    
                    for j=(obj.Solution.Trajectory.Ntrans+1):(2*obj.Solution.Trajectory.Ntrans-1)
                       pointer = j-obj.Solution.Trajectory.Ntrans+1;
                       obj.Solution.Trajectory.Body_Set(pointer)=obj.Solution.Trajectory.Body_Set(pointer).Sphere_Of_Influence();
                      SphereI = obj.Solution.Trajectory.Body_Set(pointer).SpoI;
              %         square2=(obj.Solution.Trajectory.Hyperbola(Best,pointer).Per/SphereI)^2;
              %         ceq(j)=square2-1;
                        if(bitand(obj.Solution.Trajectory.NO_ENCOUNTER,2^pointer))
                            ceq(j)=0;
                        else
                            ceq(j) =obj.Solution.Trajectory.Hyperbola(Best,pointer).Per- SphereI;
                        end
                     
                    end
                    ceqold=ceq;               
                end
                
                DeltaV=DeltaVold;

                gradient=[];
                tlast1 = tin1;
                return;

        end
        
        % Periapsis Constraints
        
        function [  con, gradient ] = Per_NLopt(tin2,run_mode)
%# Per_NLopt                :   Calculates Periapsis Constraints for Optimize_Mission

            if ~isequal(tin2,tlast1)
                    NIP=-1;
                    for j=1:obj.Solution.Trajectory.Nbody
                        if obj.Solution.Trajectory.Body_Set(j).Fixed_Point==1
                            NIP=NIP+2;
                            coslong=cos(tin1(obj.Solution.Trajectory.Nbody+NIP));
                            sinlong=sin(tin1(obj.Solution.Trajectory.Nbody+NIP));
                            coslat=cos(tin1(obj.Solution.Trajectory.Nbody+NIP+1));
                            sinlat=sin(tin1(obj.Solution.Trajectory.Nbody+NIP+1));
                            obj.Solution.Trajectory.Body_Set(j).ephem0.r(1)=obj.Min_Per(j)*coslat*coslong;
                            obj.Solution.Trajectory.Body_Set(j).ephem0.r(2)=obj.Min_Per(j)*coslat*sinlong;
                            obj.Solution.Trajectory.Body_Set(j).ephem0.r(3)=obj.Min_Per(j)*sinlat;
                        end
                    end
                    [obj.Solution, DeltaV] = obj.Solution.Compute_DeltaV(tin2(1:obj.Solution.Trajectory.Nbody));
                      DeltaVold=DeltaV;
            end
                      % Select Optimal Permutation
                 Best =  obj.Solution.Trajectory.Best;
 
                 VIOLATION = 0;
                 if (obj.Min_Per(1)>0.0)
                     ceq(1)= obj.Solution.Trajectory.dV(Best,1)/1000 - sqrt(obj.Min_Per(1));
                 end
                for j = 2:obj.Solution.Trajectory.Ntrans
                        if(bitand(obj.Solution.Trajectory.NO_ENCOUNTER,2^j))
                            ceq(j)=0;
                        else
                            ceq(j) = obj.Solution.Trajectory.Hyperbola(Best,j).Planet.radius +obj.Min_Per(j)- obj.Solution.Trajectory.Hyperbola(Best,j).Per;
              %              if ceq(j) > 0 
              %                  DeltaV = 1e50*ceq(j)/(obj.Solution.Trajectory.Hyperbola(Best,j).Planet.radius +obj.Min_Per(j));
               %                 DeltaVold=DeltaV;
                %            end
                          %  ceq(j)
                        end
                       
             %          obj.Solution.Trajectory.Body_Set(j)=obj.Solution.Trajectory.Body_Set(j).Sphere_Of_Influence();
             %          SphereI = obj.Solution.Trajectory.Body_Set(j).SpoI;
             %           square1=((obj.Solution.Trajectory.Hyperbola(Best,j).Planet.radius +obj.Min_Per(j))/SphereI)^2;
             %           square2=(obj.Solution.Trajectory.Hyperbola(Best,j).Per/SphereI)^2;
             %           ceq(j)=square1-sign(obj.Solution.Trajectory.Hyperbola(Best,j).Per)*square2;
             %          if ceq(j) >0
             %              VIOLATION=1;
             %          end
                end
                
                for j=(obj.Solution.Trajectory.Ntrans+1):(2*obj.Solution.Trajectory.Ntrans-1)
                       pointer = j-obj.Solution.Trajectory.Ntrans+1;
                       obj.Solution.Trajectory.Body_Set(pointer)=obj.Solution.Trajectory.Body_Set(pointer).Sphere_Of_Influence();
                       SphereI = obj.Solution.Trajectory.Body_Set(pointer).SpoI;
              %         square2=(obj.Solution.Trajectory.Hyperbola(Best,pointer).Per/SphereI)^2;
               %        ceq(j)=square2-1;
                         if(bitand(obj.Solution.Trajectory.NO_ENCOUNTER,2^pointer))
                            ceq(j)=0;
                         else
                            ceq(j) =obj.Solution.Trajectory.Hyperbola(Best,pointer).Per- SphereI;
                         end
                                     
                end
                
                ceqold=ceq;

           % Periapsis Constraints
           con=ceq(run_mode);
           gradient = [];
           tlast1 = tin2;
       
        end
        
        % Perihelion Constraints

       
        function [  con, gradient ] = Perhel(tin4,run_mode)
 %# Perhel                   :   Calculates Perihelion Constraints for Optimize_Mission 
 
            if ~isequal(tin4,tlast1)
                    NIP=-1;
                    for j=1:obj.Solution.Trajectory.Nbody
                        if obj.Solution.Trajectory.Body_Set(j).Fixed_Point==1
                            NIP=NIP+2;
                            coslong=cos(tin1(obj.Solution.Trajectory.Nbody+NIP));
                            sinlong=sin(tin1(obj.Solution.Trajectory.Nbody+NIP));
                            coslat=cos(tin1(obj.Solution.Trajectory.Nbody+NIP+1));
                            sinlat=sin(tin1(obj.Solution.Trajectory.Nbody+NIP+1));
                            obj.Solution.Trajectory.Body_Set(j).ephem0.r(1)=obj.Min_Per(j)*coslat*coslong;
                            obj.Solution.Trajectory.Body_Set(j).ephem0.r(2)=obj.Min_Per(j)*coslat*sinlong;
                            obj.Solution.Trajectory.Body_Set(j).ephem0.r(3)=obj.Min_Per(j)*sinlat;
                        end
                    end
                    [obj.Solution, DeltaV] = obj.Solution.Compute_DeltaV(tin4(1:obj.Solution.Trajectory.Nbody));
                      DeltaVold=DeltaV;

            end
            
            % Select Optimal Permutation
            Best =  obj.Solution.Trajectory.Best;

            for j = 1:obj.Solution.Trajectory.Ntrans
                   if (bitand(obj.Perihelia_flag,2^j))
                       PER= obj.Solution.Trajectory.Trans_Set(j).transfer_body(obj.Solution.Trajectory.perm(Best,j)).orbit.a*(1-obj.Solution.Trajectory.Trans_Set(j).transfer_body(obj.Solution.Trajectory.perm(Best,j)).orbit.e);
                       req(j) = obj.Perihelia(j+1)-PER;
                   else
                       obj.Solution.Trajectory.Trans_Set(j)=obj.Solution.Trajectory.Trans_Set(j).Calculate_Perihelion();
                       req(j) = obj.Perihelia(j+1)-obj.Solution.Trajectory.Trans_Set(j).perihelion(obj.Solution.Trajectory.perm(Best,j));
                   end
                %   if j==3
               %        obj.Solution.Trajectory.Trans_Set(j).true_anom_arr(obj.Solution.Trajectory.perm(Best,j))
              %         obj.Solution.Trajectory.Trans_Set(j).true_anom_dep(obj.Solution.Trajectory.perm(Best,j))
             %          obj.Solution.Trajectory.Trans_Set(j).perihelion(obj.Solution.Trajectory.perm(Best,j))/obj.AU
           %            obj.Perihelia(j+1)/obj.AU
            %       end
            end
                
            % Perihelion Constraints
           con=req(obj.Per_Pointer(run_mode)-1); 
          % con/obj.AU
           gradient = [];
           tlast1 = tin4;
       
        end

       
        function [cond,gradient] = Overall_Duration(tin3)
%# Overall_Duration         :   Calculates Overall Mission Duration Constraints for Optimize_Mission 
                
                cond = 0;
                for j = 2:obj.Solution.Trajectory.Nbody
                    cond=cond+tin3(j);
                end
                cond = cond-obj.Max_Duration;
             
            gradient = [];
            
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

    tplot = zeros(1,numdata);
    tt=zeros(1,numdata);
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
    
   [outstring, newpos]  = textwrap( u, D , 200);

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
    
    nplanets =EncMiss.Nbody;
    ntrans =EncMiss.Ntrans;

    % Specify Time Range

    tt=zeros(1,numdata);
    tt2=zeros(1,numdata);
    XI=zeros(numdata,3);
    XO=zeros(numdata,3);
    RI=zeros(1,numdata);
    RO=zeros(1,numdata);
    PROBIN = Body;
    PROBOU = Body;
    datei = strings(numdata);
    dateo = strings(numdata);
    
    for i=2:ntrans
        if(bitand(EncMiss.NO_ENCOUNTER,2^i))
            continue;
        end
        j=i-1;
        Best_Perm =EncMiss.perm(EncMiss.Best,i)
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
    
    D{1} = ""
    D{2} = ""
    D{3} = "   Number      Planet       Periapsis Time              Periapsis      Arr Ecc    Dep Ecc     Inclination         LOAN        AOP";
    D{4} = "                                                           km                                   degs             degs        degs";
    D{5} = ""; 
    
    for i = 1:OrbitM.Trajectory.Nbody
        Time = cspice_et2utc(OrbitM.Absolute_Times(i),'C',0);
        if i==1
            index=1;
            Data = "N/A";
        elseif(bitand(OrbitM.Trajectory.NO_ENCOUNTER,2^i))
            Data = "N/A";
        else
            index=OrbitM.Trajectory.perm(OrbitM.Trajectory.Best,i-1);
            
        end
        if i==OrbitM.Trajectory.Nbody
            index2=OrbitM.Trajectory.perm(OrbitM.Trajectory.Best,i-1);
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
    
   [outstring, newpos]  = textwrap( u, D , 200);

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
    
   [outstring, newpos]  = textwrap( u, E , 150);
    
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

    tplot = zeros(1,numdata);
    tt=zeros(1,numdata);
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
   
   MAXY1ABS= max(abs(MAXY1),abs(MINY1))
   MAXY2ABS= max(abs(MAXY2),abs(MINY2))
   MAXY3ABS= max(abs(MAXY3),abs(MINY3))
   MAXTOT=max(MAXY1ABS,MAXY2ABS)
   
 
    figure('Position', [0 0 900 900]);
    titlestr='';
    for i=1:nplanets
        titlestr=sprintf('%s %s',titlestr,PlotMiss.Body_Set(i).name);
    end
    

    for i=1:nplanets
        axis([-1.5*MAXTOT 1.5*MAXTOT -1.5*MAXTOT 1.5*MAXTOT]);
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
    for i=1:ntrans
        interval= int64(MAXFLIGHTTIME/obj.Current_Mission.Mission_Times(i+1));
        k(i)=0;
        
        datestr=cspice_et2utc(TT(i,1),'C',0);
        planetstr = PlotMiss.Body_Set(i).name;
        legendstr= sprintf('%s %s',planetstr,datestr);
        
        Best_Perm =PlotMiss.perm(PlotMiss.Best,i);
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
                p(l)=plot(-Xp2(i,l,j),Xp1(i,l,j),'or');
                XAN = [(0.5-Xp2(i,l,j)/3.0/MAXTOT-0.0001) (0.5-Xp2(i,l,j)/3.0/MAXTOT)];
                YAN = [(0.5+Xp1(i,l,j)/3.0/MAXTOT-0.0001) (0.5+Xp1(i,l,j)/3.0/MAXTOT)];

                a(l)=annotation('textarrow',XAN,YAN,'HeadStyle','none','String',PlotMiss.Body_Set(l).name,'Color','red','FontSize',9);
                hold on;
            end
            
            
            
            plot(-Y2(i,1:j),Y1(i,1:j),'Color','black');
            p(nplanets+1)=plot(-Y2(i,j),Y1(i,j),'o','Color','black');
        %    plot3(-Y2(i,1:j),Y1(i,1:j),Y3(i,1:j));
            axis([-1.5*MAXTOT 1.5*MAXTOT -1.5*MAXTOT 1.5*MAXTOT]);
       %     axis([-0.25*MAXTOT 0.25*MAXTOT -0.25*MAXTOT 0.25*MAXTOT -1e12 1e12]);
       %     axis([-1.5*MAXTOT 1.5*MAXTOT -1.5*MAXTOT 1.5*MAXTOT -1.5*MAXTOT 1.5*MAXTOT]);
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
    

