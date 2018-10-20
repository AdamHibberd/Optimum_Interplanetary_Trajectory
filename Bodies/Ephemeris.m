%# Ephemeris is a class that defines the position and velocity state of a Celestial Body
classdef Ephemeris
%#
%# PROPERTIES:
%#
%# t        : Time
%# r        : Position vector at time t (assuming cartesians)
%# R        : Radial Distance i.e. magnitude of r or norm(r)
%# v        : Velocity vector at time t (assuming cartesians)
%# V        : Speed i.e. magnitude of v or norm(v)
%#
%# METHODS:
%#
%# Ephemeris        : Class constructor can construct arrays of Ephemeris
%#

    
properties

   
    t;          %# t = Time
    r;          %# r = Position vector at time t;
    R;          %# R = Radial distance i.e. magnitude of r;
    v;          %# v = velocity vector at time t;
    V;          %# V = Speed i.e. magnitude of v;

end %# properties

methods
     
%%# Construct method for Ephemeris - Cobrtucts arrays if Necessay
function obj = Ephemeris(value)
%# Ephemeris constructs arrays of Ephemeris objects
%#
%# INPUT:
%#
%# value        : Number of Ephemeris objects to Construct
%#
%# OUTPUT:
%#
%# obj          : Array of Ephemeris objects with value number of objects
%#

    if nargin>0
        n = value;
        obj(1:n)=Ephemeris;
        for i=1:n
            obj(i).t=0;
        end
    end
end
end %# methods
    
end

