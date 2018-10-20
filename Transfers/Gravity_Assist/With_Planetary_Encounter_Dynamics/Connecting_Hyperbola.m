%# Connecting Hyperbola is a class defining 2 hyperbolic orbits which connects a Departure Velocity to an Arrival Velocity
%# assuming an impulsive velocity increment at Periapsis
classdef Connecting_Hyperbola
%#
%# PROPERTIES:
%#
%# VA       : Arrival Velocity
%# VD       : Departure Velocity 
%# alpha    : Angle between velocities VA & VD (calculated by Compute_Hyperbola)
%# alpha_thresh : Closest possible displacement of alpha from pi
%# Planet   : The Celestial Body about which the hyperbola is described
%# Probe    : The Body which is travelling along the hyperbola
%# Per      : Periapsis, Distance of closest approach by Probe to Planet's Centre (calculated by Compute_Hyperbola) 
%# beta     : Angle in Beta-plane defining plane of hyperbola (calculated by Compute_Hyperbola) 
%# DV       : Delta-V required at Periapsis to Connect VD to VA (calculated by Compute_Hyperbola) 
%# Trans1   : Transformation Matrix
%# Trans2   : Transformation Matrix
%# BTrans   : Overall Transformation Matrix
%# LOW_PER_ERROR : Periapsis below acceptable bounds
%# NO_ENCOUNTER : No Encounter dynmaics modelled for this SSO
%# NIT      : Number of Iterations Of FZERO 
%# METHODS:
%#
%# Connecting_Hyperbola     : Class Constructor initialises 2 or 1 dimensional arrays
%# Transform                : Calculates Trans1 & Trans2 based on arrival velocity vector, VA
%# Compute_Hyperbola        : Compute Hyperbola connecting Arrival Velocity VA with Departure Velocity VD
%#                          : and find Periapsis Distance, Per, DeltaV required at Periapsis and Angle beta
%# func                     : Function to be solved s.t. func=0 at solution for required trajectory
%# Calculate_Departure_Velocity     : Calculates departure velocity for zero errors required.
%#
properties
    
    VA;% Arrival Velocity
    VD;% Required Departure Velocity
    alpha;% Angle between VD & VA
    alpha_thresh = 10*pi/180; % threshold on alpha
    Planet = Body;% Central Body
    Probe = Body;% Body travelling on hyperbola wrt Central Body (Arriving)
    Probe2= Body;% Body travelling on hyperbola wrt Central Body (Departing)
    Per;% Periapsis Distance
    beta; % Angle in B-Plane, Beta
    DV;% DeltaV required at Periapsis
    Trans1; % Transformation Matrices Based on Arrival Velocity
    Trans2; % Transformation Matrices Based on Arrival Velocity
    BTrans; % Overall Transformation Matrix
    LOW_PER_ERROR; % FLAG for Low Periapsis
    NO_ENCOUNTER; % FLAG for no encounter dynamics at SSO
    NIT;		 % Number of Iterations Of FZERO 
    
end %# properties
    
methods

%%# Connecting_Hyperbola constructs 1 and 2-dim arrays of Connecting_Hyperbolas
function obj = Connecting_Hyperbola(value1,value2)
%# Connecting_Hyperbola constructs 1 and 2-dim arrays of Connecting_Hyperbolas
%#
%# INPUT:
%#
%# value1       : value1 is size of one dimensional array
%# value2       : value2 is size of 2nd dimension of 2 dimensional array if present
%#
%# OUTPUT:
%#
%# obj          : Connecting_Hyperbolas with arrays of VA & VD initialised.
%#

if nargin>1
        obj(1:value1,1:value2)=Connecting_Hyperbola;
        for i=1:value1
            for j=1:value2
                obj(i,j).VA = zeros(1,3);
                obj(i,j).VD = zeros(1,3);
            end
        end
    elseif nargin == 1
        obj(1:value1)=Connecting_Hyperbola;
        for i=1:value1
            obj(i).VA = zeros(1,3);
            obj(i).VD = zeros(1,3);
        end
    else
        obj.VA = zeros(1,3);
        obj.VD = zeros(1,3);
end
    

end %# Connecting_Hyperbola

%%# Transform calculated Transformation Matrices from Ecliptic to B-Plane base on Arrival Velocity
function obj = Transform(obj)
%# Transform calculated Transformation Matrices
%#
%# INPUT:
%#
%# obj          : Object in question with uninitialsied Transformation Matrices
%#
%# OUTPUT:
%#
%# obj          : Object in question with initialsied Transformation Matrices
%#
    % Projection of VA onto ecliptic
    Vij = sqrt(obj.VA(2)^2+obj.VA(1)^2);

    % Magnitude of VA (speed)
    Vmag = norm(obj.VA);

    % Compute angles

    cosdel = obj.VA(2)/Vij;
    sindel = obj.VA(1)/Vij;
    cosgam = obj.VA(3)/Vmag;
    singam = Vij/Vmag;

    % Transformation Matrices
  %  obj.Trans1 = [ cosdel, sindel, 0; -sindel, cosdel, 0; 0, 0, 1 ];
    obj.Trans1 = [ cosdel, -sindel, 0; sindel, cosdel, 0; 0, 0, 1 ];
%    obj.Trans2 = [ 1, 0, 0 ; 0, cosgam, singam ; 0, -singam, cosgam];
    obj.Trans2 = [ 1, 0, 0 ; 0, cosgam, -singam ; 0, singam, cosgam];

end   %# Transform

%%# Compute_Hyperbola connecting Arrival Velocity VA with Departure Velocity VD
%%# and find Perapsis Distance , DeltaV required at Periapsis and Angle beta
 function obj = Compute_Hyperbola(obj)
%# Compute_Hyperbola connecting Arrival Velocity VA with Departure Velocity VD
%#
%# INPUT:
%#
%# obj          : Object in question with given VA, VD
%#
%# OUTPUT:
%#
%# obj          : Object in question with calculated alpha, beta, Per and DV
%#


    % Determine angle alpha between Departure Velocity and Arrival
     % Velocity

     cosalpha = dot(obj.VA, obj.VD ) / norm(obj.VA) / norm(obj.VD) ;
     obj.alpha = acos(cosalpha);
     VAV = sqrt(norm(obj.VA)*norm(obj.VD));

     
     if (isempty(obj.Planet.mu)|| obj.Planet.mu==0)
         obj.Per=[];
         obj.DV = norm(obj.VD-obj.VA);
         obj.NO_ENCOUNTER = 1;
     elseif (abs(obj.alpha)>pi-obj.alpha_thresh)
         obj.Per = obj.Planet.mu*(pi-abs(obj.alpha))^2/8/VAV^2;
         obj.LOW_PER_ERROR = 1;
         obj.NO_ENCOUNTER = 0;
         obj.DV =  abs(1-obj.Per/obj.Planet.radius)*1e50 ;
     else
         obj.LOW_PER_ERROR = 0;
         obj.NO_ENCOUNTER = 0;
  

    % vtemp1 = -cross(obj.VA,obj.VD);
    % vtemp2 = [ obj.VA(1)*obj.VA(3), obj.VA(3)*obj.VA(2), -(obj.VA(1)^2 + obj.VA(2)^2)];

     % Determine angle beta

 %    cosbeta = dot(vtemp1,vtemp2)/norm(obj.VA)^2/norm(obj.VD)/sqrt(obj.VA(1)^2 + obj.VA(2)^2)/sin(obj.alpha);
 %    obj.beta = 2*pi-acos(cosbeta)

  %  obj.beta =  acos(cosbeta);
  
      sinbeta = (obj.VA(2)*obj.VD(1) - obj.VA(1)*obj.VD(2))/norm(obj.VD)/sqrt(obj.VA(1)^2 + obj.VA(2)^2)/sin(obj.alpha); 
      obj.beta =  asin(sinbeta);
       % obj.beta*180/pi
      % Guess Angle between periapsis and Arrival Asymptote

       X0 = (pi + obj.alpha )/2;
    
        obj.NIT = 0;
        X = FZERO(X0,1e-50, 0, 200);
        if (obj.NIT>=200)
            X0 = pi + obj.alpha - X;
            X = FZERO(X0,1e-50,1,400);
            X = pi + obj.alpha - X;
        end

    %   X
    %   Vmin2 = min( norm(obj.VD)^2,norm(obj.VA)^2);
    %   P0 = - obj.Planet.mu/VAV * ( 1 + 1/cos(X0));
   %     P0=1.4720e9;
  %  [P fval exitflag output ] = fzero(@func2,P0,options);

      Per = -obj.Planet.mu/norm(obj.VA)^2 * ( 1 + 1/cos(X));
  %      Per = P;

     obj.DV = sqrt(2*obj.Planet.mu/Per + norm(obj.VD)^2) - sqrt(2*obj.Planet.mu/Per + norm(obj.VA)^2) ;

     obj.Per = Per;
     end
  
      % Vdep = Calculate_Departure_Velocity(obj)

     return;
     
%# Employs Newton Iteration to converge to solution for angle theta (output x)
function x = FZERO(x0, xtol, mode,  Maxit)


	reduce = 1.0;

	funcold = 0;

	VD2 = norm(obj.VD)^2;
	VA2 = norm(obj.VA)^2;
	cal = cos(obj.alpha);
	sal = sin(obj.alpha);
	
	x = x0;
	for i=obj.NIT:Maxit
	
		if (mode == 0)
		
			func = (VA2 + VD2*cal)*cos(x) + (VD2 - VA2)*sal*sin(x)*cos(x) + (VD2 - VA2)*cos(x)^2 * cal + VD2*sal*sin(x);
			func = func / obj.Planet.mu;
		%	func = 1 / obj.Planet.mu * (((VD2 - VA2)*cos(x) + VD2)*cos(obj.alpha - x) + VA2*cos(x));
			dfuncbydx = 1/obj.Planet.mu*((VD2 - VA2)*sin(obj.alpha - 2 * x) + VD2 * sin(obj.alpha - x) - VA2 *sin(x));
		
		else
			
			func = ((VA2 - VD2)*cos(obj.alpha - x) + VD2) *cos(x) + VA2*cos(obj.alpha-x);
			dfuncbydx = (VA2 - VD2)*sin(obj.alpha - 2 * x) + VA2 * sin(obj.alpha - x) - VD2 *sin(x);
		
        end
		if ((funcold * func) < 0.0)
            
            reduce = reduce*0.1;

		else
	
			reduce = 1.0;
        end
		
		funcold = func;
		xold = x;

		dx = reduce*func / dfuncbydx;
		x = x - dx;		

		if (abs(x-xold) < xtol) 
            break;
        end
    end
	obj.NIT = i;

	if (x >  2*pi)
        
        x = x -  2*pi*int(x / 2 / pi);
    elseif (x < -2 * pi)
	
		x = abs(x);
		x = pi - x + pi*int(x / pi);
    end
	return;
end


%%# Calculate_Departure_Velocity based on solution calculated by Compute_Hyperbola (should equal VD)
function Vdep = Calculate_Departure_Velocity(obj)
%# Calculate_Departure_Velocity based on solution calculated by Compute_Hyperbola (should equal VD)
%#
%# INPUT:
%#
%# obj          : Connecting_Hyperbola Object with Computed Hyperbolas
%#
%# OUTPUT:
%#
%# Vdep         : Departure Velocity should equal VD
%#

% Compute Arrival Speed

     Vmag = norm(obj.VA);

     % First Compute Velocity at Periapsis

     Periapsis = obj.Per;
     Beta = obj.beta;
     DV = obj.DV;

     VPer = sqrt(2*obj.Planet.mu/Periapsis + Vmag^2);

     % Calculate Arrival Eccentricity

     e = Periapsis*Vmag^2/obj.Planet.mu + 1 ;

     % Calculate Angle between Arrival Velocity and Periapsis

     costheta = -1/e ;

     sintheta = sqrt(1-costheta^2);

     % Calculate Departing Eccentricity

     edash = Periapsis * (VPer+DV)^2 / obj.Planet.mu - 1;

     % Calculate Angle Between Departure Velocity and Periapsis

     costhetadash = -1/edash ;

     sinthetadash = sqrt(1-costhetadash^2);

     % Work out the Departure Velocity Try in the B-Plane First

     Vdep2 = (VPer+DV)^2 - 2*obj.Planet.mu/Periapsis;

     Vdepmag = sqrt( (VPer+DV)^2 - 2*obj.Planet.mu/Periapsis);
     Vdep =  [ 0 , -Vdepmag*(sintheta*costhetadash+costheta*sinthetadash) , -Vdepmag*(sintheta*sinthetadash-costhetadash*costheta) ];

     % Now Rotate to the Ecliptic

     T = [ cos(Beta), sin(Beta), 0 ; -sin(Beta), cos(Beta), 0 ; 0 , 0 , 1 ];

     Vtemp =  -transpose(obj.Trans1) * transpose(obj.Trans2) * transpose(T) * transpose(Vdep);
     Vdep = transpose(Vtemp);

     return;
     
end %#Calculate_Departure_Velocity

 end %# Compute_Hyperbola
 
function obj = Orbits_From_Hyperbolas(obj)

    % Initialise Gravitational Masses
    obj.Probe.orbit.GM = obj.Planet.mu;
    obj.Probe2.orbit.GM = obj.Planet.mu;
    
    % Arrival Speed
    VAS = norm(obj.VA);
    % Departure Speed
    VDS = norm(obj.VD);
    
    % Arrival Energy
    ENA = VAS*VAS/2;
    % Departure Energy
    END = VDS*VDS/2;
    
    % Arrival Semi-major Axis
    obj.Probe.orbit.a = -obj.Probe.orbit.GM/2/ENA;
    obj.Probe.orbit.arec = 1/obj.Probe.orbit.a;
    obj.Probe.orbit.arec = -ENA*2/obj.Probe.orbit.GM;
    % Departure Semi-major Axis
    obj.Probe2.orbit.a = -obj.Probe2.orbit.GM/2/END;
    obj.Probe2.orbit.arec = 1/obj.Probe2.orbit.a;
     obj.Probe2.orbit.arec = -END*2/obj.Probe.orbit.GM;
     
    % Arrival Eccentricity
    obj.Probe.orbit.e = obj.Per*VAS*VAS/obj.Probe.orbit.GM + 1;
    % Departure Eccentricity
    obj.Probe2.orbit.e = obj.Per*VDS*VDS/obj.Probe2.orbit.GM + 1;
    
    % Arrival Semi-latus Rectum or Parameter
    obj.Probe.orbit.p = obj.Probe.orbit.a*(1-obj.Probe.orbit.e^2);
    % Departure Semi-latus Rectum or Parameter
    obj.Probe2.orbit.p = obj.Probe2.orbit.a*(1-obj.Probe2.orbit.e^2); 
    
    % Arrival Angular Momentum
    HA = sqrt(obj.Probe.orbit.p*obj.Probe.orbit.GM);
    % Departure Angular Momentum
    HD = sqrt(obj.Probe2.orbit.p*obj.Probe2.orbit.GM);
    
    % Unit Angular Momentum Vector
    uH = cross(obj.VA,obj.VD)/VAS/VDS/sin(obj.alpha);
    
    % Arrival Angular Momenum Vector
    HAV = HA * uH;
    % Departure Angular Momentum Vector
    HDV = HD * uH;
         
    % Longitude of Ascending Node
    L = atan2( uH(1) , -uH(2) );
    obj.Probe.orbit.loan = L;
    obj.Probe2.orbit.loan = L;
    
    % Inclination
    I = atan2( sqrt( uH(1)^2 + uH(2)^2 ) , uH(3) );
    obj.Probe.orbit.I = I;
    obj.Probe2.orbit.I = I;
    
    % Laplace-Range-Lenz vector for Arrival EvA
    EvA = cross( obj.VA , HAV ) + obj.Probe.orbit.GM/VAS * obj.VA;
    EA = norm(EvA);
    % Laplace-Range-Lenz vector for Departure EvD
    EvD = cross( obj.VD , HDV ) - obj.Probe2.orbit.GM/VDS * obj.VD;
    ED = norm(EvD);

    % Arrival  Argument of Perigee wA needs to be calculated
    sinwA = EvA(3) / EA / sin(I);
    coswA = ( EvA(1) + EA * sinwA * cos(I) * sin(L) ) / EA / cos(L);
    obj.Probe.orbit.aop = atan2( sinwA , coswA );
    
    % Departure Argument of Perigee wD needs to be calculated
    sinwD = EvD(3) / ED / sin(I);
    coswD = ( EvD(1) + ED * sinwD * cos(I) * sin(L) ) / ED / cos(L);
    obj.Probe2.orbit.aop = atan2( sinwD , coswD );

    % Initialise True Anomalies to equal zero
    obj.Probe.orbit.ta0 = 0;
    obj.Probe2.orbit.ta0 = 0;
    
    % Inialise Epochs to Zero
    obj.Probe.orbit.epoch = 0;
    obj.Probe2.orbit.epoch = 0;
    
end
     
 
end %# methods
end

