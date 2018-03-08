%# Mission is a class for determining Interplanetary Trajectories based on given Bodies and given time inputs
classdef Mission
%#
%# PROPERTIES:
%#
%#        
%#        Trajectory;           : Trajectory Of Mission
%#        Mission_Times;        : Array of relative Encounter Times (first element is Absolute Time)
%#        Absolute_Times;       : Array of Absolute Encounter Times
%#        Spice_Min_Times;      : Array of Minimum Times for SPICE Calculation
%#        Spice_Max_Times;      : Array of Maximum Times for SPICE Calculation
%#        Out_Of_Spice_Bounds =0; : Flag to indicate SPICE was unable to be used for a body i
%#        TotaldV;              : Total DeltaV of Mission
%#        FlybyRendez=0;        : Flyby or Rendezvous Flag =0 (FLYBY) = 1 (RENDEZVOUS)
%#        wayflag=1;            : Default to Prograde Only
%#
%# METHODS:
%#
%#
%# Mission                  : Class Constructor
%# Set_Mission_Times        : Function to Calculate Mission Times from input array of Absolute Times
%# Set_Absolute_Times       : Function to Calculate Absolute Times from input array of Mission Times
%# Compute_DeltaV           : Function to Update Trajectory Based on input Mission Times and output DeltaV
%#
    properties
    
        Trajectory;             % Trajectory Of Mission
        Mission_Times;          % Array of relative Encounter Times
        Absolute_Times;         % Array of Absolute Encounter Times
        Spice_Min_Times;        % Array of Minimum Times for SPICE Calculation
        Spice_Max_Times;        % Array of Maximum Times for SPICE Calculation
        Out_Of_Spice_Bounds =0; % Flag to indicate SPICE was unable to be used for a body i
        TotaldV;                % Total DeltaV of Mission
        FlybyRendez=0;          % Flyby or Rendezvous Flag =0 (FLYBY) = 1 (RENDEZVOUS)
        wayflag=1;              % Default to Prograde Only
 
    end
    
    methods
        
        function obj = Mission( bodies , times , MinSpice, MaxSpice)
%# Mission                  : Class Constructor

            % Initialise SPICE Min & MAx Times
            obj.Spice_Min_Times = MinSpice
            obj.Spice_Max_Times = MaxSpice
            
            % Initialise Positions and Osculating Orbits of Various Bodies
 
            for i=1:size(bodies,2)
                if times(i)<obj.Spice_Min_Times(i)
                    mode(i)=1;
                    obj.Out_Of_Spice_Bounds=bitor(obj.Out_Of_Spice_Bounds,2^i,'uint32');
                    bodies(i)=bodies(i).compute_ephem_at_t(obj.Spice_Min_Times(i),2,1e-4);
                    bodies(i)=bodies(i).calculate_orbit_from_ephem(obj.Spice_Min_Times(i));
                    bodies(i).ephem0=bodies(i).ephemt;
                elseif times(i)>obj.Spice_Max_Times(i)
                    mode(i)=1;
                    obj.Out_Of_Spice_Bounds=bitor(obj.Out_Of_Spice_Bounds,2^i,'uint32');
                    bodies(i)=bodies(i).compute_ephem_at_t(obj.Spice_Max_Times(i),2,1e-4);
                    bodies(i)=bodies(i).calculate_orbit_from_ephem(obj.Spice_Max_Times(i));
                    bodies(i).ephem0=bodies(i).ephemt;
                
                else
                    if bitand(obj.Out_Of_Spice_Bounds,2^i)
                        obj.Out_Of_Spice_Bounds=obj.Out_Of_Spice_Bounds-2^i;
                    end
                    mode(i)=2;
                    bodies(i)=bodies(i).compute_ephem_at_t(times(i),mode(i),1e-4);
                    bodies(i)=bodies(i).calculate_orbit_from_ephem(times(i));
                end

            end
            
            % Create an object of Type Nbody_Trajectory_With_Encounters
            obj.Trajectory = Nbody_Trajectory_With_Encounters(bodies);
        

            obj.Mission_Times = zeros(1,obj.Trajectory.Nbody);
            obj.Absolute_Times = zeros(1,obj.Trajectory.Nbody);
            
            obj = obj.Set_Mission_Times( times );
            
            [obj DeltaV gradient] = obj.Compute_DeltaV( obj.Mission_Times );
            
        end
        
        function obj = Set_Mission_Times( obj, times )
%# Set_Mission_Times        : Function to Calculate Mission Times from input array of Absolute Times

            % Function to Initialise Mission_Times 
            % Input times is array of absolute times
            % Output has first item as launch time and
            % remaining times as cruise durations
            
            obj.Mission_Times(1) = times(1);
            
            for i=2:obj.Trajectory.Nbody
                obj.Mission_Times(i)=times(i)-times(i-1);
            end
            
        end
        
        function obj = Set_Absolute_Times( obj, times )
 %# Set_Absolute_Times       : Function to Calculate Absolute Times from input array of Mission Times
 
            % Reverse of Set_Mission_Times
            
            obj.Absolute_Times(1)=times(1);
            
            for i=2:obj.Trajectory.Nbody
                obj.Absolute_Times(i)=times(i)+obj.Absolute_Times(i-1);
            end
        
        end
        
        function [obj, DeltaV, gradient] = Compute_DeltaV( obj, times )
%# Compute_DeltaV           : Function to Update Trajectory Based on input Mission Times and output DeltaV

            obj = obj.Set_Absolute_Times(times);
            for i=1:obj.Trajectory.Nbody
               
                if obj.Absolute_Times(i)<obj.Spice_Min_Times(i)
                    mode(i)=1;
                    obj.Out_Of_Spice_Bounds=bitor(obj.Out_Of_Spice_Bounds,2^i,'uint32');
                    obj.Trajectory.Body_Set(i)=obj.Trajectory.Body_Set(i).compute_ephem_at_t(obj.Spice_Min_Times(i),2,1e-4);
                    obj.Trajectory.Body_Set(i)=obj.Trajectory.Body_Set(i).calculate_orbit_from_ephem(obj.Spice_Min_Times(i));
                    obj.Trajectory.Body_Set(i).ephem0=obj.Trajectory.Body_Set(i).ephemt;
                elseif obj.Absolute_Times(i)>obj.Spice_Max_Times(i)
                    mode(i)=1;
                    obj.Out_Of_Spice_Bounds=bitor(obj.Out_Of_Spice_Bounds,2^i,'uint32');
                    obj.Trajectory.Body_Set(i)=obj.Trajectory.Body_Set(i).compute_ephem_at_t(obj.Spice_Max_Times(i),2,1e-4);
                    obj.Trajectory.Body_Set(i)=obj.Trajectory.Body_Set(i).calculate_orbit_from_ephem(obj.Spice_Max_Times(i));
                    obj.Trajectory.Body_Set(i).ephem0=obj.Trajectory.Body_Set(i).ephemt;
                else
                    mode(i)=2;
                end
            end
            obj.Trajectory = obj.Trajectory.Compute_Total_Deltav( obj.Absolute_Times, mode, 1e-4, 1000,obj.wayflag,obj.FlybyRendez);
            DeltaV = obj.Trajectory.BestDeltaV;
            gradient=[];
            obj.TotaldV = DeltaV;
        end
            
            
    end
end

