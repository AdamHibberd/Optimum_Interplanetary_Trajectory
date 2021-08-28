%# Definition for Orbit Class - includes standard orbital parameters
classdef Orbit
%#
%# PROPERTIES:
%#
%# epoch        : Time of orbital parameters including reference true anomaly, ta0
%# ta0          : Reference true anomaly at epoch, ta0
%# ta           : Current true anomaly
%# p            : Semi-latus Rectum (otherwise known as parameter)
%# a            : Semi-major axis
%# arec         : Reciprocal of semi-major axis, 1/a 
%# e            : Eccentricity
%# aop          : Argument of Perigee (radians)
%# loan         : Longitude of Ascending Node (radians)
%# I            : Inclination (radians)
%# GM           : Gravitational Mass of Attracting Centre (default = Sun)
%# TP           : Time Period of Orbit
%# Trans_PtoE   : Perifocal to Ecliptic Transformation Matrix
%#
%# METHODS:
%#
%# get.Trans_PtoE   : Computes Trans_PtoE based on orbital parameters
%#    
        
properties
    
    epoch;              %# epoch = time of orbital parameters
    ta0;                %# ta0 = true anomaly at epoch
    ta;               %# ta = true anomaly
    p;                  %# p = parameter or semi-latus rectum
    a;                  %# a = semi-major axis
    arec;               %# arec = 1/a
    e;                  %# e = eccentricity
    aop;                %# aop = argument of perigee
    loan;               %# loan = longitude of ascending node
    I;                  %# I = inclination
    GM = 1.32712440018e20; %# GM = gravitational mass of attracting body (Assume Sun)
    TP;                 %# TP = Orbital Time period for elliptical orbits
    Trans_PtoE;         %# Trans_PtoE = Perifocal to Ecliptic Transformation Matrix            

end %# properties
    
methods
    
%%# Function get.Trans_PtoE computes the Perifocal to Ecliptic Transformation Matrix based on Orbital Parameters
function val = get.Trans_PtoE(obj)
%# get.Trans_PtoE computes the Perifocal to Ecliptic Transformation Matrix 
%#
%# INPUT:
%# None
%#
%# OUTPUT:
%#
%# val          : Trans_PtoE = Perifocal to Ecliptic Transformation Matrix
%#
    Trans_OtoE(1,1) = cos(obj.loan);
    Trans_OtoE(1,2) = -sin(obj.loan)*cos(obj.I);
    Trans_OtoE(1,3) = sin(obj.I)*sin(obj.loan);
    Trans_OtoE(2,1) = sin(obj.loan);
    Trans_OtoE(2,2) = cos(obj.loan)*cos(obj.I);
    Trans_OtoE(2,3) = -sin(obj.I)*cos(obj.loan);
    Trans_OtoE(3,1) = 0;
    Trans_OtoE(3,2) = sin(obj.I);
    Trans_OtoE(3,3) = cos(obj.I);
    Trans_PtoE = zeros(3,3);
    Trans_PtoO(1,1) = cos(obj.aop);
    Trans_PtoO(1,2) = -sin(obj.aop);
    Trans_PtoO(1,3) = 0;
    Trans_PtoO(2,1) = sin(obj.aop);
    Trans_PtoO(2,2) = cos(obj.aop);
    Trans_PtoO(2,3) = 0;
    Trans_PtoO(3,1) = 0;
    Trans_PtoO(3,2) = 0;
    Trans_PtoO(3,3) = 1;

    val = Trans_OtoE * Trans_PtoO;
    
end %# get.Trans_PtoE


end %# methods
    
end

