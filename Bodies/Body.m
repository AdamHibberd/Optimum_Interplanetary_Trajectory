%# Body is a class that defines a Celestial Body
classdef Body
%#
%# PROPERTIES:
%#
%# name     : string identifier for Celestial Body, eg 'Earth'
%# ID       : integer defining Celestial Body SPICE ID
%# time     : Current time (secs)
%# mu       : Gravitational Mass of Body (m3/s2)
%# GM       : Gravitational Mass of Attracting Centre (Sun) (m3/s2)
%# radius   : Equatorial radius of Body (m)
%# orbit    : Orbit of Celestial Body (Orbit object)
%# ephem0   : Ephemeris of Celestial Body at Epoch (normally t=0) (Ephemeris Object)
%# ephemt   : Ephemeris of Celestial Body at t = time (Ephemeris Object)
%# SPICEFLAG: = 0 Indicates SPICE successfully Calculated Ephemeris
%#          : = -1 Indicates SPICE failure and Alternative Orbit Calculation
%#           : Prefered
%# state    : state of Body as furnished by SPICE -   state(1:3) = position (km)
%#                                                    state(4:6) = velocity (km/s)
%# 
%# METHODS:
%#
%# Body         : Constructor can construct arrays of Body
%# orbittoephem : Calculates ephemeris ephem0 using Orbit parameters based upon the true anomaly
%# Cz           : Special Function of Auxiliary Variable z (Universal Variable Formulation)
%#                - see page 208 of Fundamentals of Astrodynamics, Bate, Mueller, White
%# Sz           : Special Function of Auxiliary Variable z (Universal Variable Formulation)
%#                - see page 208 of Fundamentals of Astrodynamics, Bate, Mueller, White
%# dCdz         : Derivative of Cz wrt z
%# dSdz         : Derivative of Sz wrt z
%# compute_ephem_at_t : Calculates ephemeris ephemt at time t using Orbit parameters : method used depends upon mode
%# calculate_orbit_from_ephem : Calculates orbit from ephemeris ephemt
%# Sphere_Of_Influence : Calculates Laplace Sphere of Influence for bosy/Sun Combination
%# Calculate_True_Anomaly : Calculates True Anomaly
    
properties
    
        name;       %# string identifier for Celestial Body, eg 'Earth'
        ID;         %# integer defining Celestial Body SPICE ID
        time;       %# Current time (secs)
        mu;         %# Gravitational Mass of Body (m3/s2)
        GM;         %# Gravitational Mass of Attracting Centre (Sun) (m3/s2)
        radius;     %# Equatorial radius of Body (m)
        orbit = Orbit; %# Orbit of Celestial Body (Orbit object)
        ephem0 = Ephemeris; %# Ephemeris of Celestial Body at Epoch (normally t=0)
        ephemt = Ephemeris; %# Ephemeris of Celestial Body at t = time
        state;      %# state of Body as furnished by SPICE
        SpoI;       %# Sphere of Influence as suggested by Laplace
        Fixed_Point = 0; %# =0 Means Body is in orbit, =1 Means Intermediate Point, Otherwise a Fixed Point
        
end %# properties

methods
        
%%# Constructor Method for Body  - Constructs arrays if necessary
function obj = Body(value)
%# Body Constructs arrays of Body objects
%#
%# INPUT:
%#
%# value        : Number of Body objects to Construct
%#
%# OUTPUT:
%#
%# obj          : Array of Body objects with value number of objects
%#
    if nargin>0
        n = value;
        obj(1:n)=Body;
        for i=1:n
            obj(i).time=0;
            obj(i).orbit = Orbit;
            obj(i).ephem0 = Ephemeris;
            obj(i).ephemt = Ephemeris;
        end
    else
        obj.time=0;
        obj.orbit = Orbit;
        obj.ephem0 = Ephemeris;
        obj.ephemt = Ephemeris;
    end
end %# Body
                    
%%# Calculate Ephemeris from Orbit
function obj = orbittoephem(obj,t)
%# orbittoephem Calculates ephemeris ephem0 using Orbit parameters based upon the true anomaly
%#
%# INPUT:
%#
%# obj          : Current Body object in question
%# t            : Epoch for ephemeris ephem0.t
%#
%# OUTPUT:
%#
%# obj          : Body object with Ephemeris ephem0 calculated from orbit and true anomaly ta0 
%#            
    % Time for ephemeris data
    obj.ephem0.t = t;
    
    % Radial distance
    obj.ephem0.R = obj.orbit.p / ( 1 + obj.orbit.e * cos(obj.orbit.ta0) );
            
    % Position components in perifocal system
    obj.ephem0.r(1,1) = obj.ephem0.R * cos(obj.orbit.ta0) ;
    obj.ephem0.r(2,1) = obj.ephem0.R * sin(obj.orbit.ta0) ;
    obj.ephem0.r(3,1) = 0;
            
    % Position components in ecliptic system
    obj.ephem0.r = obj.orbit.Trans_PtoE * obj.ephem0.r ;
            
    % Radial velocity
    radialv =  sqrt(obj.orbit.GM/obj.orbit.p) * obj.orbit.e * sin(obj.orbit.ta0);
            
    % Horizontal velocity in perifocal system
    horizv = sqrt(obj.orbit.GM/obj.orbit.p) * (1 + obj.orbit.e * cos(obj.orbit.ta0) );
            
    % Velocity components in perifocal system
    obj.ephem0.v(1,1) = -sqrt(obj.orbit.GM/obj.orbit.p)*sin(obj.orbit.ta0);
    obj.ephem0.v(2,1) =  sqrt(obj.orbit.GM/obj.orbit.p)*(obj.orbit.e + cos(obj.orbit.ta0) );
    obj.ephem0.v(3,1) = 0;
    obj.ephem0.V = sqrt(radialv^2 + horizv^2);
            
    % Velocity components in ecliptic system
    obj.ephem0.v = obj.orbit.Trans_PtoE * obj.ephem0.v ;
   
end %# orbittoephem
    
%%# Function Cz using universal variable formulation
function val = Cz(obj,z)
%# Cz Special Function of Auxiliary Variable z
%#
%# INPUT:
%#
%# obj          : Current Body object in question
%# z            : Auxiliary variable z
%#
%# OUTPUT:
%#
%# val          : value of Cz function at z 
%#        

    if abs(z) < 0.001
        val = 1/2 - z/4/3/2 + z^2/6/5/4/3/2 - z^3/8/7/6/5/4/3/2;
    elseif z <= -0.001
        val = ( 1 - cosh(sqrt(-z)) ) / z ;
    else
        val = ( 1 - cos(sqrt(z)) ) / z ;
    end
end %# Cz
        
%%# Function Sz using universal variable formulation        
function val = Sz(obj,z)
%# Sz Special Function of Auxiliary Variable z
%#
%# INPUT:
%#
%# obj          : Current Body object in question
%# z            : Auxiliary variable z
%#
%# OUTPUT:
%#
%# val          : value of Sz function at z 
%#  
    if abs(z) < 0.001
        val = 1/3/2 - z/5/4/3/2 + z^2/7/6/5/4/3/2 - z^3/9/8/7/6/5/4/3/2;
    elseif z <= -0.001
        val = ( sinh(sqrt(-z)) - sqrt(-z) ) / sqrt(-z)^3 ;
    else
        val = ( sqrt(z) - sin(sqrt(z)) ) / sqrt(z)^3 ;
    end
end %# Sz
                
%%# Function dCdz using universal variable formulation
function val = dCdz(obj,z)
%# dCdz Derivative of Cz wrt z
%#
%# INPUT:
%#
%# obj          : Current Body object in question
%# z            : Auxiliary variable z
%#
%# OUTPUT:
%#
%# val          : value of Cz gradient function at z 
%#
    if abs(z) < 0.001
        val = -1/4/3/2 + 2*z/6/5/4/3/2 -3*z^2/8/7/6/5/4/3/2 + 4*z^3/10/9/8/7/6/5/4/3/2;
    else
        val = 1/2/z*(1 - z*obj.Sz(z) - 2*obj.Cz(z)) ;
    end
end %# dCdz
        
%%# Function dSdz using universal variable formulation        
function val = dSdz(obj,z)
%# dSdz Derivative of Sz wrt z
%#
%# INPUT:
%#
%# obj          : Current Body object in question
%# z            : Auxiliary variable z
%#
%# OUTPUT:
%#
%# val          : value of Sz gradient function at z 
%#
    if abs(z) < 0.001
        val = -1/5/4/3/2 + 2*z/7/6/5/4/3/2 - 3*z^2/9/8/7/6/5/4/3/2 + 4*z^3/11/10/9/8/7/6/5/4/3/2;
    else
        val = 1/2/z*(obj.Cz(z) - 3*obj.Sz(z)) ;
    end
end %# dSdz
        
%%# Compute ephem at time t
function obj = compute_ephem_at_t(obj, t, mode, thresh)
%# compute_ephem_at_t Calculates ephemeris ephemt at time t using Orbit parameters : method used depends upon mode
%#
%# INPUT:
%#
%# obj          : Current Body object in question
%# t            : Time at which ephemeris ephemt.t will be calculated
%# mode         : mode = 1 use Universal Formulation from Fixed Orbital Parameters
%#              : mode = 2 Use SPICE
%# thresh       : time tolerance for numerical iteration required for mode = 1
%#
%# OUTPUT:
%#
%# obj          : Body object with Ephemeris ephemt calculated from orbit at time t 
%#            
    obj.time = t;

    if (obj.Fixed_Point>0)
        obj.ephemt = obj.ephem0;
        obj.ephemt.t = t;
        obj.ephemt.R = norm(obj.ephemt.r);
        obj.ephemt.V = norm(obj.ephemt.v);
        return;
    else
        obj.ephemt.t = t;
    end
        
    
    if ( mode==2 )
        state = cspice_spkezr(obj.ID,t,'ECLIPJ2000','NONE','SUN');

        obj.ephemt.r(1:3) = state(1:3)*1e+3;
        obj.ephemt.v(1:3) = state(4:6)*1e+3;
        obj.ephemt.R = norm(obj.ephemt.r);
        obj.ephemt.V = norm(obj.ephemt.v);
        return;
                
    else
                
        % Overall change in time
        Dt = obj.time - obj.ephem0.t;
        if (abs(Dt) <= thresh)
            obj.ephemt = obj.ephem0 ;
            return;
        end
        % Correct for Orbital Time Period if necessary
            
        % Dt = Dt - obj.orbit.TP*fix(Dt/obj.orbit.TP);
            
        % r0.v0 = initial dot product of r0 and v0
        rdotv0 = dot(obj.ephem0.r, obj.ephem0.v) ;
            
        % Initial Guess for x (universal variable)
        if obj.orbit.e <= 1
            xn = sqrt(obj.orbit.GM)*obj.orbit.arec*Dt;
        else
            xn = sign(Dt)*sqrt(-1/obj.orbit.arec)*log(-2*obj.orbit.arec*obj.orbit.GM*Dt/(rdotv0 + sign(Dt)*sqrt(-obj.orbit.GM/obj.orbit.arec)*(1-obj.orbit.arec*obj.ephem0.R)));
        end
            
        zn =    xn^2*obj.orbit.arec;
        tn =    rdotv0*xn^2*obj.Cz(zn)/obj.orbit.GM + (1 - obj.orbit.arec*obj.ephem0.R)*xn^3*obj.Sz(zn)/sqrt(obj.orbit.GM) + obj.ephem0.R*xn/sqrt(obj.orbit.GM);
        dtdxn = rdotv0*xn*(1 - zn*obj.Sz(zn))/obj.orbit.GM + xn^2*obj.Cz(zn)/sqrt(obj.orbit.GM) + obj.ephem0.R * (1 - zn*obj.Cz(zn))/sqrt(obj.orbit.GM);
            
        % Do Newton iteration
        i=0;
        while abs(tn - Dt) > thresh
            i=i+1;
            if i>1000
                break;
            end
            xn = xn + (Dt - tn) / dtdxn ;
            zn =    xn^2*obj.orbit.arec;            
            tn =    rdotv0*xn^2*obj.Cz(zn)/obj.orbit.GM + (1 - obj.orbit.arec*obj.ephem0.R)*xn^3*obj.Sz(zn)/sqrt(obj.orbit.GM) + obj.ephem0.R*xn/sqrt(obj.orbit.GM);
            dtdxn = rdotv0*xn*(1 - zn*obj.Sz(zn))/obj.orbit.GM + xn^2*obj.Cz(zn)/sqrt(obj.orbit.GM) + obj.ephem0.R * (1 - zn*obj.Cz(zn))/sqrt(obj.orbit.GM);               
        end

        % Compute f & g for position vector at time Dt

        f = 1 - xn^2*obj.Cz(zn)/obj.ephem0.R;
        g = Dt - xn^3*obj.Sz(zn)/sqrt(obj.orbit.GM);

        obj.ephemt.r = f * obj.ephem0.r + g * obj.ephem0.v ;
        obj.ephemt.R = norm(obj.ephemt.r);

        % Compute fdot & gdot for velocity vector at time Dt

        gdot = 1 - xn^2*obj.Cz(zn)/obj.ephemt.R;
        fdot = sqrt(obj.orbit.GM)/obj.ephemt.R/obj.ephem0.R*xn*(zn*obj.Sz(zn) - 1 );

        obj.ephemt.v = fdot * obj.ephem0.r + gdot * obj.ephem0.v ;
        obj.ephemt.V = norm(obj.ephemt.v);

        % Compute true anomaly

        obj.orbit.ta = atan2( dot(obj.ephemt.r,obj.ephemt.v)/obj.ephemt.R*sqrt(obj.orbit.p/obj.orbit.GM), obj.orbit.p/obj.ephemt.R - 1 );

        % Correct for Orbital Time Period for Elliptical Orbits
    end

end %# compute_ephem_at_t
        
%%# calculate_orbit_from_ephem : Calculates orbit from ephemeris ephemt        
function obj = calculate_orbit_from_ephem( obj , time)
%# Calculates osculating orbital parameters from cartesians
%#
%# INPUT:
%#
%# obj          : Current Body object in question
%# t            : Epoch for ephemeris ephemt.t
%#
%# OUTPUT:
%#
%# obj          : Body object with orbit calculated from ephemeris data 
%#    
        % Gravitational Mass
        mu = obj.orbit.GM;

        % Epoch for orbit

        obj.orbit.epoch= time;

        % Radial Distance

        x = obj.ephemt.r;
        v = obj.ephemt.v;
        R = norm(x);

        % Speed

        V = norm(v);

        % Angular Momentum Vector h

        h = cross( x , v );

        % Magnitude of angular momentum H

        H = norm(h); 

        % Longitude of Ascending Node

        L = atan2( h(1) , -h(2) );

        obj.orbit.loan = L;

        % Inclination

        I = atan2( sqrt( h(1)^2 + h(2)^2 ) , h(3) );
        obj.orbit.I = I;

        % Laplace-Range-Lenz vector Ev

        Ev = cross( v , h ) - mu/R * x;

        E = norm(Ev);

        % Argument of Perigee w needs to be calculated
        
        sinw = Ev(3) / E / sin(I);

        cosw = ( Ev(1) + E * sinw * cos(I) * sin(L) ) / E / cos(L);

        % w = atan2( Ev(3) * cos(L) , ( sin(I)* Ev(1) + Ev(3) * cos(I) * sin(L) ) );
        obj.orbit.aop = atan2(sinw, cosw);

        % Energy per unit mass gives semi-major axis a

        Energy = V * V / 2 - mu/R;

        if abs(Energy) > 0 

            obj.orbit.a = -mu / 2 / Energy;
        end;

        obj.orbit.arec = -Energy*2/mu;

        % Eccentricity, e 

        e = sqrt( 1 + 2 * Energy * ( H / mu )^2 );
        obj.orbit.e= e;

        % Parameter, p (semi=latus rectum)

        obj.orbit.p = H^2/mu;

        % Time Period
        if ( e < 1)

            obj.orbit.TP = 2*pi*sqrt(obj.orbit.a^3/mu);
        end

            % True Anomaly

        rxd = x(1) * cos(L) + x(2) * sin(L);
        ryd = x(2) * cos(L) - x(1) * sin(L);
        rydd = ryd * cos(I) + x(3) * sin(I);

    %    sinta = ( rydd * cosw - rxd * sinw ) / R;
    %    costa = ( rydd * sinw - rxd * cosw ) / R;

    %    true_anom = atan2( sinta , costa );
        
 %       Y = rydd * Ev(1)* sin(I) + rydd * Ev(3) * cos(I)*sin(L) - rxd *Ev(3)*cos(L);
 %       X = rydd * Ev(3)* cos(L) - rxd * Ev(3) * cos(I)*sin(L)  - rxd * Ev(1)*sin(I);
        
 %       true_anom = atan2( Y , X );
 %       sinta = ( rydd * cos(w) - rxd * sin(w) ) / R;
 %       costa = ( rydd * sin(w) + rxd * cos(w) ) / R;

 %       true_anom = atan2( sinta , costa );
 %       obj.orbit.ta = true_anom;
        obj = obj.Calculate_True_Anomaly();

end %# calculate_orbit_from_ephem

% Calculate Sphere of Influence

function obj = Sphere_Of_Influence(obj)

    obj.SpoI = obj.orbit.a*(obj.mu/obj.orbit.GM)^(0.4);

end

% Calculate the True Anomaly


function obj = Calculate_True_Anomaly(obj)
    
    % Gravitational Mass
    mu = obj.orbit.GM;
         
    % Radial Distance

    x = obj.ephemt.r;
    v = obj.ephemt.v;
    R = norm(x);
    
    % Angular Momentum Vector h

    h = cross( x , v );

    % Magnitude of angular momentum H

    H = norm(h); 

    % Parameter, p (semi=latus rectum)

    p = H^2/mu;
    
    Radial_Velocity= dot(x,v)/R;
    
    
    if (R==p)&&(Radial_Velocity>0)
        true_anom=pi/2;
    elseif (R==p)&&(Radial_Velocity<0)
        true_anom=-pi/2;
    else
        true_anom=atan2(R*Radial_Velocity*sqrt(p/mu),(p-R));
    end
   
    obj.orbit.ta = true_anom;
    
end
    
    


end %# methods



end 

