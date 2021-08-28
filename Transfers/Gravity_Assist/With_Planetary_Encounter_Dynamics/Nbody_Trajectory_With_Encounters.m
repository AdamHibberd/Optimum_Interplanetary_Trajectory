%# Nbody_Trajectory_With_Encounters is a subclass of Nbody_Trajectory
%# It includes a 2D array of hyperbolas connecting for each Body a Departure
%# Velocity with an Arrival Velocity
classdef Nbody_Trajectory_With_Encounters < Nbody_Trajectory
%#
%# PROPERTIES:
%#
%# Hyperbola    : 2D Array of Connecting_Hyperbola 
%#
%# METHODS:
%#
%# Nbody_Trajectory_With_Encounters     : Class Constructor initialises 2 dimensional arrays of Hyperbola
%# Compute_Total_Deltav                 : Computes Overall combined DeltaV at each Periapsis
%#

properties

    DELV;           %# Array of DeltaV's
    HYPERB;         %# Array of Hyperbolas
    Hyperbola;      %# Connecting Hyperbolas
    DIVERGING;      %# Interplanetary Trajectory Has unacceptably low Periapsis here
    NO_ENCOUNTER;   %# Interplanetary Trajectory has no encounter dynamics for this SSO

end %# properties

methods

%%# Nbody_Trajectory_With_Encounters is class constructor
function obj = Nbody_Trajectory_With_Encounters(bodies)
%%# Nbody_Trajectory_With_Encounters is class constructor
%#
%# INPUT:
%#
%# bodies       : List of Celestial Bodies to be connected (each of Body class)
%#
%# OUTPUT:
%#
%# obj          : The Nbody_Trajectory in question, Initialised
%#

        obj = obj@Nbody_Trajectory(bodies);
        obj.Hyperbola = Connecting_Hyperbola(obj.NP,obj.Ntrans);
        obj.HYPERB = Connecting_Hyperbola(4,obj.Ntrans);

end %# Nbody_Trajectory_With_Encounters

%%# Compute_Total_Deltav computes Overall DeltaV for different times for a
%%# multi-planet Trajectory with Planetary Encounters
function obj = Compute_Total_Deltav(obj ,t, mode, thresh, maxit, wayflag,flybyrendez,home_periapsis,target_periapsis)
%# Compute_Total_Deltav computes Overall DeltaV for different times for a multi-planet Trajectory with Planetary Encounters
%# INPUT:
%#
%# obj          : Current Nbody_Trajectory in question
%# t            : Array of times t(1)= Launch Time but t(2)etc are traj durations
%# array of modes : Mode of operation for calculating Ephemeris (see Body class method compute_ephem_at_t)
%# thresh       : Threshold for calculating Ephemeris (see Body class method compute_ephem_at_t)
%# maxit        : Maximum number of iterations for Calculate_transfer
%# flybyrendez  : Flyby or rendezvous with destination planet
%#
%# OUTPUT:
%#
%# obj          : The Nbody_Trajectory_With_Encounters in question, with all the Delta V variables calculated
%#    

   obj = Compute_Total_Deltav@Nbody_Trajectory(obj ,t, mode, thresh, maxit, wayflag, flybyrendez,home_periapsis,target_periapsis);

   obj.DeltaV(1:obj.NP)=0.0;

   if (wayflag==0)

        for j=2:obj.Ntrans

                obj.HYPERB(1,j).VA(:)=obj.VA(:,j,1);
                obj.HYPERB(1,j).VD(:)=obj.VD(:,j,1);
                obj.HYPERB(1,j).Planet = obj.Body_Set(j);
           %     obj.HYPERB(1,j)=obj.HYPERB(1,j).Transform();
                obj.HYPERB(1,j)=obj.HYPERB(1,j).Compute_Hyperbola();

                
                obj.HYPERB(2,j).VA(:)=obj.VA(:,j,2);
                obj.HYPERB(2,j).VD(:)=obj.VD(:,j,1);
                obj.HYPERB(2,j).Planet = obj.Body_Set(j);
          %      obj.HYPERB(2,j)=obj.HYPERB(2,j).Transform();
                obj.HYPERB(2,j)=obj.HYPERB(2,j).Compute_Hyperbola();
            
                
                obj.HYPERB(3,j).VA(:)=obj.VA(:,j,1);
                obj.HYPERB(3,j).VD(:)=obj.VD(:,j,2);
                obj.HYPERB(3,j).Planet = obj.Body_Set(j);
         %       obj.HYPERB(3,j)=obj.HYPERB(3,j).Transform();
                obj.HYPERB(3,j)=obj.HYPERB(3,j).Compute_Hyperbola();
             
                
                obj.HYPERB(4,j).VA(:)=obj.VA(:,j,2);
                obj.HYPERB(4,j).VD(:)=obj.VD(:,j,2);
                obj.HYPERB(4,j).Planet = obj.Body_Set(j);
         %       obj.HYPERB(4,j)=obj.HYPERB(4,j).Transform();
                obj.HYPERB(4,j)=obj.HYPERB(4,j).Compute_Hyperbola();
        end
   
        for i=1:obj.NP
            for j=1:obj.Nbody
                if (j<=obj.Ntrans)&&(obj.Trans_Set(j).ni(obj.perm(i,j)))>=maxit
                    obj.DeltaV(i)=inf;
                    continue;
                end
                  if j==1
                      if home_periapsis>0
                        VelPer=sqrt(2*obj.Body_Set(j).mu/home_periapsis+norm(obj.VD(:,j,obj.perm(i,j)))^2)-sqrt(obj.Body_Set(j).mu/home_periapsis);
                      else
                        VelPer= norm(obj.VD(:,j,obj.perm(i,j)));
                      end
                        obj.dV(i,j) = VelPer;
                  elseif j==obj.Nbody
                      if flybyrendez>0
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
                     if (obj.perm(i,j)==1) && (obj.perm(i,j-1)==1)
                        obj.Hyperbola(i,j) = obj.HYPERB(1,j);
                     elseif (obj.perm(i,j)==1) && (obj.perm(i,j-1)==2)
                        obj.Hyperbola(i,j) = obj.HYPERB(2,j);
                     elseif (obj.perm(i,j)==2) && (obj.perm(i,j-1)==1)
                        obj.Hyperbola(i,j) = obj.HYPERB(3,j);
                     elseif (obj.perm(i,j)==2) && (obj.perm(i,j-1)==2)
                        obj.Hyperbola(i,j) = obj.HYPERB(4,j);
                     end
                     obj.dV(i,j) = abs(obj.Hyperbola(i,j).DV);
                  end
            obj.DeltaV(i)=obj.DeltaV(i)+obj.dV(i,j);
            end
        end
        [obj.BestDeltaV, obj.Best] = min(obj.DeltaV,[],2);

   
   else
        i = obj.Best;

        for j=1:obj.Nbody
               if j==1
                      if home_periapsis>0
                        VelPer=sqrt(2*obj.Body_Set(j).mu/home_periapsis+norm(obj.VD(:,j,obj.perm(i,j)))^2)-sqrt(obj.Body_Set(j).mu/home_periapsis);
                      else
                        VelPer= norm(obj.VD(:,j,obj.perm(i,j)));
                      end
                      obj.dV(i,j) = VelPer;
               elseif j==obj.Nbody
                      if flybyrendez>0
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
                        obj.Hyperbola(i,j).VA(:) = obj.VA(:,j,obj.perm(i,j-1));
                        obj.Hyperbola(i,j).VD(:) = obj.VD(:,j,obj.perm(i,j));
                        obj.Hyperbola(i,j).Planet = obj.Body_Set(j);
                   %     obj.Hyperbola(i,j)=obj.Hyperbola(i,j).Transform();
                        obj.Hyperbola(i,j)=obj.Hyperbola(i,j).Compute_Hyperbola();
                        obj.dV(i,j) = abs(obj.Hyperbola(i,j).DV);
               end
               obj.DeltaV(i)=obj.DeltaV(i)+obj.dV(i,j);
        end

        obj.BestDeltaV = obj.DeltaV(i);
 
   end
   
   obj.DIVERGING=0;
   obj.NO_ENCOUNTER=0;
    if obj.Nbody>2
         for j=2:obj.Ntrans
                if (obj.Hyperbola(obj.Best,j).LOW_PER_ERROR ==1)
                        obj.DIVERGING=obj.DIVERGING + 2^j;
                end
                if (obj.Hyperbola(obj.Best,j).NO_ENCOUNTER ==1)
                        obj.NO_ENCOUNTER=obj.NO_ENCOUNTER + 2^j;
                end
         end
    end

 
end  %# Compute_Total_Deltav
end %# methods
end

