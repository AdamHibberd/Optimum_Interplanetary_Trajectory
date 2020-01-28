%# Nbody_Trajectory is a class defining a set of trajectories connecting 3 or more Celestial Bodies at specified times
classdef Nbody_Trajectory
%#
%# PROPERTIES:
%#
%# Name     : Name of Nbody_Trajectory (Mission Name)
%# Nbody    : Number of Celestial Bodies to be connected by transfers
%# Ntrans   : Number of transfers i.e. Nbody - 1 
%# NP       : Number of permutations for connecting transfers
%#          : i.e. if each transfer has a long and short ways there are 
%#          : NP = (2 raised to power of Ntrans) possibilities
%# perm     : Array of all possible permutations (NP elements)
%# Body_Set : List of Nbody Celestial Bodies to be Connected by transfers
%# Trans_Set: List of Ntrans Transfer_orbits, each connecting 2 adjacent Celestial Bodies 
%#          : in Body_Set and each with long way and short way
%# DeltaV   : Array of Total Connecting DeltaVs, one for each of the NP permutations
%# BestDeltaV: Mimimum Value element of the DeltaV array i.e. Best Total DeltaV
%# Best     : Best Permutation from DeltaV array
%# deltaV   : Array of delta v vectors, 3 dimensional array i.e. (3 * Ntrans * 2)
%#          : i.e. there are Ntrans * 2 velocity changes required (long way and short way) 
%#          : 1st is Departure velocity required at Home planet (First Planet in list)
%#          : The remaining are changes required at intermediate planets
%#          : Assume no velocity change for final planet, so not included in array
%# dV       : Array of Magnitudes of deltaV's 2 Dimensional array i.e. (NP * Ntrans) 
%# VD     : Array of Departure Velocities, 3 Dimensional array (3 * Nbody * 2) 
%# VA     : Array of Arrival Velocities,3 Dimensional array (3 * Nbody * 2)  
%# 
%# METHODS:
%#
%# Nbody_Trajectory         : Class Constructor accepts as input a list of Celestial Bodies to be connected 
%# Initialise_Transfers     : Initialises transfer arrays and variables
%# Calculate_Nbody_ephem    : Computes the Ephemeris for all of the Celestial Bodies           : Special Function of Auxiliary Variable z (Universal Variable Formulation)
%# Calculate_Nbody_transfers: Calculates the transfer orbits connecting up all the Bodies     : Special Function of Auxiliary Variable z (Universal Variable Formulation)
%# Compute_Total_Deltav     : Calculates delta v arrays
%#  
 
properties  
        
    Name;           %# Name of Nbody_Trajectory Object
    Nbody;          %# Number of Celestial Bodies to visit including first and last body
    Ntrans;         %# Number of transfers
    NP;             %# Possible Permutations
    perm;           %# Array of All Permutations
    Body_Set;       %# List of Bodies to Visit
    Trans_Set;      %# List of Intermediate Transfers
    DeltaV;         %# Total Connecting DeltaVs
    BestDeltaV;     %# Best Total DeltaV
    Best;           %# Best Permutation
    deltaV;         %# Array of deltav vectors
    dV;             %# Array of Magnitudes of deltaV
    VD;           %# Array of Departure Velocities
    VA;           %# Array of Arrival Velocities
    
end %# properties

methods

%%# Nbody_Trajectory is the class constructor method 
function obj = Nbody_Trajectory(bodies)
%# Nbody_Trajectory is the class constructor
%#
%# INPUT:
%#
%# bodies       : List of Celestial Bodies to be connected (each of Body class)
%#
%# OUTPUT:
%#
%# obj          : The Nbody_Trajectory in question, Initialised
%#
   if nargin>0
   % Initialise set of NBodies for connection

   obj.Nbody = size(bodies,2);
   obj.Ntrans = obj.Nbody-1;
   obj.Body_Set = Body(obj.Nbody);
   
   for i=1:obj.Nbody
       obj.Body_Set(i) = bodies(i);
   end

   % Initialise set of Transfer Orbits

   obj.Trans_Set = Transfer_orbit(obj.Ntrans);


   obj = Initialise_Transfers(obj);

   % Initialise DeltaVs            
   obj.deltaV=zeros(3,obj.Ntrans,2);
   obj.dV=zeros(obj.NP,obj.Nbody);
   obj.DeltaV=zeros(1,obj.NP);
   obj.BestDeltaV = 0;
   obj.VD = zeros(3,obj.Nbody,2);
   obj.VA = zeros(3,obj.Nbody,2);
   else
       return;
   end
end %# Nbody_Trajectory

%%# Initialise_Transfers initialises permutation information and transfer orbits
function obj = Initialise_Transfers(obj)
%# Initialise_Transfers initialises permutation information and transfer orbits
%#
%# INPUT:
%#
%# obj          : Current Nbody_Trajectory in question
%#
%# OUTPUT:
%#
%# obj          : The Nbody_Trajectory in question, Initialised
%#
   % Initialise Number of Permutations
   obj.NP = 2^obj.Ntrans;    
   for bytes = 0:obj.NP-1
         for i =1:obj.Ntrans
             obj.perm(bytes+1,i)= 1 + bitand(bytes,2^(i-1))/2^(i-1);  
         end
   end    

    % Initialise set of Transfer Orbits

   for i=1:obj.Ntrans
       obj.Trans_Set(i).bodyd = obj.Body_Set(i);
       obj.Trans_Set(i).bodya = obj.Body_Set(i+1);
       for j=1:2
            obj.Trans_Set(i).transfer_body(j).ephem0 = obj.Body_Set(i).ephemt;
            obj.Trans_Set(i).transfer_body(j).ephemt = obj.Body_Set(i+1).ephemt;
       end
   end

end %# Initialise_Transfers

%%# Calculate_Nbody_ephem computes list of Ephemeris for each body of the Body_Set
function obj = Calculate_Nbody_ephem(obj,mode,thresh)
%# Calculate_Nbody_ephem computes list of Ephemeris for each body of the Body_Set
%#
%# INPUT:
%#
%# obj          : Current Nbody_Trajectory in question
%# mode         : Array of Modes of operation for calculating Ephemeris (see Body class method compute_ephem_at_t)
%# thresh       : Threshold for calculating Ephemeris (see Body class method compute_ephem_at_t)
%#
%# OUTPUT:
%#
%# obj          : The Nbody_Trajectory in question, with all the Ephemeris of the Bodies calculated
%#

    for i = 1:obj.Nbody
        obj.Body_Set(i)=obj.Body_Set(i).compute_ephem_at_t(obj.Body_Set(i).ephemt.t,mode(i),thresh);
    end
    
end %# Calculate_Nbody_ephem

%%# Calculate_Nbody_transfers computes list of Transfer Orbits Connecting each of the Bodies successively
function obj = Calculate_Nbody_transfers(obj, thresh, itmax, wayflag)
%# Calculate_Nbody_transfers computes list of Transfer Orbits Connecting each of the Bodies successively
%#
%# INPUT:
%#
%# obj          : Current Nbody_Trajectory in question
%# thresh       : Threshold for calculating Transfer Orbits (see Transfer_orbit class method Calculate_transfer)
%# itmax        : Maximum number of iterations for Calculate_transfer
%#
%# OUTPUT:
%#
%# obj          : The Nbody_Trajectory in question, with all the Transfer_orbits calculated
%#

% Compute list of Transfer Orbits

    for i=1:obj.Ntrans
        td = obj.Trans_Set(i).bodyd.ephemt.t;
        tar = obj.Trans_Set(i).bodya.ephemt.t;
        obj.Trans_Set(i)=obj.Trans_Set(i).Calculate_transfer(td,tar,thresh,itmax,wayflag);
    end
end %# Calculate_Nbody_transfers

%%# Compute_Total_Deltav computes Overall DeltaV for different times for a multi-planet Trajectory 
function obj = Compute_Total_Deltav(obj ,t, mode, thresh, maxit,wayflag, flybyrendez,home_periapsis,target_periapsis)
%# Compute_Total_Deltav computes Overall DeltaV for different times for a multi-planet Trajectory
%#
%# INPUT:
%#
%# obj          : Current Nbody_Trajectory in question
%# t            : Array of times t(1)= Launch Time but t(2)etc are traj durations
%# mode         : Array of Mode of operation for calculating Ephemeris (see Body class method compute_ephem_at_t)
%# thresh       : Threshold for calculating Ephemeris (see Body class method compute_ephem_at_t)
%# maxit        : Maximum number of iterations for Calculate_transfer
%# flybyrendez  : Flyby or rendezvous with destination planet
%#
%# OUTPUT:
%#
%# obj          : The Nbody_Trajectory in question, with all the Delta V variables calculated
%#

    for i = 1:obj.Nbody
        obj.Body_Set(i).ephemt.t = t(i);
     %   if i>1
     %       obj.Body_Set(i).ephemt.t = obj.Body_Set(i-1).ephemt.t+t(i);
     %   end
    end

    obj = obj.Calculate_Nbody_ephem(mode,thresh);  

    obj = obj.Initialise_Transfers();

    obj = obj.Calculate_Nbody_transfers(thresh, maxit, wayflag);

    for j = 1:obj.Nbody
        for k = 1:2
            if j==1
                % Start with Departure Planet
                obj.VA(:,j,k)=0;
                obj.VD(:,j,k) = obj.Trans_Set(j).transfer_body(k).ephem0.v - obj.Body_Set(j).ephemt.v;
                obj.deltaV(:,j,k) = obj.VD(:,j,k);
                continue;
            end
            if j==obj.Nbody
                % End With Final Planet
                if (flybyrendez==0)
                    obj.VD(:,j,k) = obj.Trans_Set(j-1).transfer_body(k).ephemt.v - obj.Body_Set(j).ephemt.v;
                else
                    obj.VD(:,j,k) = 0;
                end
                obj.VA(:,j,k) = obj.Trans_Set(j-1).transfer_body(k).ephemt.v - obj.Body_Set(j).ephemt.v;
                continue;
            end
            obj.VD(:,j,k) = obj.Trans_Set(j).transfer_body(k).ephem0.v - obj.Body_Set(j).ephemt.v;
            obj.VA(:,j,k) = obj.Trans_Set(j-1).transfer_body(k).ephemt.v - obj.Body_Set(j).ephemt.v;
            obj.deltaV(:,j,k) = obj.VD(:,j,k) - obj.VA(:,j,k) ; 
        end
    end

    obj.DeltaV=zeros(1,obj.NP);

    for i = 1:obj.NP
        for j = 1:obj.Nbody
            if j==1
                 if home_periapsis>0
                      VelPer=sqrt(2*obj.Body_Set(j).mu/home_periapsis+norm(obj.VD(:,j,obj.perm(i,j)))^2)-sqrt(obj.Body_Set(j).mu/home_periapsis);
                  else
                      VelPer= norm(obj.VD(:,j,obj.perm(i,j)));
                  end
                  obj.dV(i,j) = VelPer;
             elseif j==obj.Nbody
                if(flybyrendez>0)
                    if target_periapsis>0
                        V2=sqrt(2*obj.Body_Set(j).mu/target_periapsis);
                        V1=sqrt(norm(obj.VA(:,j,obj.perm(i,j-1)))^2+V2^2);
                        VelTar=abs(V2-V1);
                    else
                        VelTar=norm(obj.VA(:,j,obj.perm(i,j-1)));
                    end
                else
                    VelTar = 0;
                end
                obj.dV(i,j)=VelTar;
            else
                obj.dV(i,j)=norm(obj.VD(:,j,obj.perm(i,j))-obj.VA(:,j,obj.perm(i,j-1)));
            end
            obj.DeltaV(i)=obj.DeltaV(i)+obj.dV(i,j);
        end
    end

[obj.BestDeltaV obj.Best] = min(obj.DeltaV,[],2);     

end %# Compute_Total_Deltav       

end %# methods

end

