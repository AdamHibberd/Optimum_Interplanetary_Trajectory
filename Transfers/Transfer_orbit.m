%# Transfer_Orbit connects 2 Bodies at specified times with 2 transfer obits, long way and short way 
classdef Transfer_orbit
%#
%# PROPERTIES:
%#
%# ni       : 2 element array corresponding to number of iterations for convergence for both ways
%# nmax     : Maximum number of iterations allowed before giving up
%# bodyd    : Departure Celestial Body - a Body object
%# bodya    : Arrival Celestial Body - a Body object
%# transfer_body : 2 Transfer Bodies (2 element array of Body objects) - One for short way, other for long way
%# ephemd   : Departure Ephemeris (2 element array of Ephemeris objects) - One for short way, other for long way
%# ephema   : Arrival Ephemeris (2 element array of Ephemeris objects) - One for short way, other for long way
%# td       : Time of Departure from Departure Body
%# tar      : Time of Arrival at Arrival Body
%# dVd      : Delta-V at Departure - for short way and long way
%# dVa      : Delta-V at Arrival - for long way and short way
%#
%# METHODS:
%#
%# Transfer_orbit       : Class constructor constructs arrays of Transfer_orbits
%# Calculate_transfer   : Calculate 2 Transfer Orbits from Departure Body to Arrival Body
%#                        Uses Universal Variables Formulation as described in Fundamentals
%#                        of Astrodynamics, Bate, Mueller, White
%# FZERO                : Function Calculates the value of z for which y is zero
%# deltaVd              : Calculate Departure DeltaVs
%# deltaVa              : Calculate Arrival DeltaVs
%#

    % Transfer Orbit from Body 1 to Body 2
    
properties

        ni= zeros(1,2);     %#  2 element array corresponding to number of iterations for convergence for both ways
        nmax;               %# Maximum number of iterations     
        bodyd;              %# Departure Planet/Body
        bodya;              %# Arrival Planet/Body
        transfer_body;      %# Transfer Body - One for short way and the other for long way
        ephemd;             %# Departure ephemeris - Long way and short way
        ephema;             %# Arrival ephemeris - Long way and Short way
        true_anom_dep;      %# True anomaly at Departure - Long way and Short Way
        true_anom_arr;      %# True anomaly at Arrival - Long way and Short Way
        perihelion;         %# Perihelia of Transfer ARC - Long Way and Short way
        td;                 %# Departure Time
        tar;                %# Arrival Time
        dVd = zeros(1,2);   %# Delta-V at Departure - for long way and short way
        dVa = zeros(1,2);   %# Delta-V at Arrival - for long way and short way


end %# properties
    
methods
    
%%# Transfer_orbit Class Constructor Method -  Construct arrays if Necessay
function obj = Transfer_orbit(value)
%# Transfer_orbit Constructs arrays of Transfer_orbit objects
%#
%# INPUT:
%#
%# value        : Number of Transfer_orbit objects to Construct
%#
%# OUTPUT:
%#
%# obj          : Array of Transfer_orbit objects with value number of objects
%#
    if nargin>0
        n = value;
        obj(1:n)=Transfer_orbit;
        for i=1:n
            obj(i).td=0;
            obj(i).ephemd = Ephemeris(2);
            obj(i).ephema = Ephemeris(2);
            obj(i).transfer_body = Body(2);
            obj(i).bodya = Body;
            obj(i).bodyd = Body;
         end
    else
        obj.ephemd = Ephemeris(2);
        obj.ephema = Ephemeris(2);
        obj.transfer_body = Body(2); 
        obj.bodya = Body;
        obj.bodyd = Body;
    end

end %# Transfer_orbit

%%# Calculate Transfer Orbit from Departure Body to Body 2
%%# Uses Universal Variables Formulation as described in Fundamentals
%%# of Astrodynamics, Bate, Mueller, White
function obj = Calculate_transfer(obj,td,tar,thresh,itmax,wayflag)
%# Calculate_transfer Calculates 2 transfer orbits (long way & short way) between 2 Bodies
%#
%# INPUT:
%#
%# obj          : Current Transfer_orbit object in question
%# td           : Departure Time from bodyd (Departure Body)
%# tar          : Arrival Time at bodya (Arrival Body)
%# thresh       : Threshold for iteration process
%# itmax        : Maximum number of iteration cycles before giving up
%# wayflag      : =0 calculate short way and longway
%#              : =1 calculate prograde only
%#
%# OUTPUT:
%#
%# obj          : Transfer_orbit object with departure and arrival ephemeris and orbit calculated
%#            

    obj.nmax = itmax;

    % Set Departure and Arrival times amd Positions
    obj.td = td;
    obj.tar = tar;
    
    for i1 =1:2
        obj.ephemd(i1) = obj.bodyd.ephemt;
        obj.ephema(i1) = obj.bodya.ephemt;
    end
    
    % Calculate unit angular momentum vector
    he=cross(obj.bodyd.ephemt.r,obj.bodyd.ephemt.v);
    he=he/norm(he);
    
    % Calculate semi-major axis
    sma = (obj.bodya.ephemt.R+obj.bodyd.ephemt.R)/2.0;


    % Departure radial distance
    rd = obj.bodyd.ephemt.R;

    % Arrival radial distance
    ra = obj.bodya.ephemt.R;

    % Sun's Gravitational Mass
    GM = obj.bodya.orbit.GM;
    
%    sma = (((tar-td)/pi)^2*GM)^(1/3);

    % Start to iterate to Solution based on j=1 (one way) then 
    %  try j=2 (opposite way)

    for j = 1:2
                 % Change in true anomaly
        cosdta = dot(obj.bodya.ephemt.r,obj.bodyd.ephemt.r)/ra/rd;

        dta = acos(cosdta);

        if j==1
            dta =2*pi-dta;
        end
              
        % Calculate constant A
        A = sqrt(rd*ra)*sin(dta)/sqrt(1-cos(dta));

        if (A >= 0 && dta <= pi)
            z0=0.0;
            [znmin ,dzn ]= FZERO(z0, 1e-1, 1000);
            zn = znmin+dzn;
        else

            % Guess inital value of z
            zn = (3/2*pi)^2;
            znmin=0.0;
        end
        speed_dep= sqrt(GM*(2*sma-rd)/sma/rd);
        speed_arr= sqrt(GM*(2*sma-ra)/sma/ra);
        dta1=100*speed_dep/rd;
        dta2=100*speed_arr/ra;
        
        if (dta < pi + dta1+dta2) && (dta > pi -dta1-dta2)
            pi_flag = 1;
            if speed_dep<speed_arr
                obj.ephemd(j).v = -speed_dep*cross(obj.ephemd(j).r,he)/rd ;
                obj.ephema(j).v = -speed_arr*cross(obj.ephema(j).r,he)/ra ;
            elseif speed_dep>speed_arr
                obj.ephemd(j).v = speed_dep*cross(obj.ephemd(j).r,he)/rd ;
                obj.ephema(j).v = speed_arr*cross(obj.ephema(j).r,he)/ra ;
            else
                obj.ephemd(j).v = speed_dep*cross(obj.ephemd(j).r,he)/rd ;
                obj.ephema(j).v = speed_arr*cross(obj.ephema(j).r,he)/ra ;
            end
            if j==2
                obj.ephemd(j).v = -obj.ephemd(j).v;
                obj.ephema(j).v = -obj.ephema(j).v;
            end
        else
            pi_flag=0;
        end
        % Calculate Special Functions of zn and their gradient wrt
        % zn
        Sn = obj.bodya.Sz(zn);
        Cn = obj.bodya.Cz(zn);
        dSdzn = obj.bodya.dSdz(zn);
        dCdzn = obj.bodya.dCdz(zn);

        % Calculate special universal variables for initial guess
        yn = ra + rd - A * ( 1 -zn*Sn) / sqrt(Cn);
        xn = sqrt( yn / Cn );
        

        % Calculate corresponding Time-of-Flight for initial guess
        tn = (xn^3 * Sn + A * sqrt(yn))/ sqrt(GM);

        % Calculate Gradient of tn wrt zn
        dtdzn = (xn^3 * (dSdzn - 3*Sn*dCdzn/2/Cn) + A/8*(3*Sn*sqrt(yn)/Cn + A/xn))/sqrt(GM);

        % Calculate initial step dzn
        dzn = (tar-td-tn)/dtdzn;

        % Inialise Error in ToF wrt Required ToF
        deltatnew=tar-td-tn;
        deltatold=deltatnew;

        % Iniialise Variables for Newton interation
        i=0;
        factor=1;
        errflag=0;
        newtry=0;
        exitflag =0;
        limitznflag=0;
        limitynflag=0;
        limittnflag=0;
        
        if ~pi_flag
            % Do Newton Iteration on zn
            while not(exitflag)

               try 
                   exitflag = (abs(tn-tar+td)<thresh);

               catch
                   exitflag=0;
               end

                i=i+1;

                % Break out if not converging
                if (i>itmax)||(exitflag)
                    break;
                end

                % Don't allow factor to get to small - try a different
                % guess for dzn
      %          if (abs(factor*dzn)<thresh/abs(dtdzn))
                 if (factor<1e-10)
                        errflag=1;
                        factor=1;
                 end
                if (limitznflag==0)&&(limitynflag==0)&&(limittnflag==0)
                    factor=1;
                end
                % Normal change in zn
                if (errflag==0)

                    dzn = factor*deltatold/dtdzn;

                % Else if Not Converging try a totally different guess
                else
                    zn=0;
     %               RAND=rand()-0.5;
                     RAND = rand();
                     if znmin==0
                         dzn= sign(RAND-0.5)*0.04*(RAND-0.5)^2*pi^2;
                     else
                         Range = (2*pi)^2-znmin;
                         dzn = Range*RAND+znmin;
                     end
     %               dzn= sign(RAND-0.5)*(0.1*(RAND-0.5)*pi)^2; 
                    errflag=0;
                    newtry=1;

                end

                % Iterate to solution
                zn=zn+dzn;

                % Don't allow zn to get to large or small
              %  for ii=1:itmax

              %  end

                % Calculate Special Functions of zn and their gradient wrt
                % zn
                Sn = obj.bodya.Sz(zn);
                Cn = obj.bodya.Cz(zn);
                dSdzn = obj.bodya.dSdz(zn);
                dCdzn = obj.bodya.dCdz(zn);


                % Calculate Universal Variable yn based on zn 
                yn = ra + rd - A * ( 1 -zn*Sn) / sqrt(Cn);

                % Don't allow yn to become -ve
                if ( yn < 0 )
                    zn = zn - dzn;
                    factor=factor*0.1;
                    limitynflag=1;
                    continue; 
                else
                    limitynflag=0;
                end

                % Calculate Universal Variable xn based on yn and zn
                xn = sqrt( yn / Cn );

                % Calculate ToF 
                tn = (xn^3 * Sn + A * sqrt(yn))/ sqrt(GM);

                % Calculate Error in ToF wrt required 

                deltatnew=tar-td-tn;


                % Determine Change in ToF wrt previous iteration
                deltatime=deltatnew-deltatold;
    %            for ii=1:itmax/4
                     if( zn > (2*pi)^2 )
                            zn=zn-dzn;
                            factor=factor*0.1;
                            limitznflag=1;
                            continue;

                     else
                         limitznflag=0;
       %                   break;
                     end
         %           break;

                if abs(deltatnew)<thresh
                    break;
                end
      %          end
       %         if ii>1
        %            continue;
         %       end
                % If ToF has reached a local minimum try a new guess
                % for dzn
                if ((abs(deltatime)<thresh/100000)) &&newtry==0

     %           if ((abs(deltatime)<thresh/100000)&&(abs(deltatime)>0.0)) &&newtry==0
    %                i
    %                deltatime
    %                deltatnew
     %               zn
      %              yn
                    errflag=1;
                    continue;
                else
                    errflag=0;
                end      

                % Don't allow the new guess to be further away from the
                % solution than the old guess
              %  if (abs(deltatold)-abs(deltatnew))<0&newtry==0
              %      '3'
              %      rateflag=1; 
              %      zn= zn - dzn;
              %     factor=factor*0.1;
              %     continue;
              %  end


                % Don't go backwards in time
                if (tn < 0)

                    zn = zn - dzn;
                    factor=factor*0.1;
                    limittnflag=1;
                    continue;
                else
                    limittnflag=0;
                end




                % Calculate Gradient of ToF wrt zn
                dtdzn = (xn^3 * (dSdzn - 3*Sn*dCdzn/2/Cn) + A/8*(3*Sn*sqrt(yn)/Cn + A/xn))/sqrt(GM);

                % Don't allow Gradient of ToF to become -ve
         %       if (dtdzn<0&newtry==0)
         %           '5'
         %           zn=zn-dzn;
         %           factor=factor*0.1;
         %            dtdzn=dtdznold;
         %       end

                % Normal Operation of iteration
                newtry=0;
                deltatold=deltatnew;
           end


           % Update Number of iterations
           obj.ni(j)=i;

           % Compute f & g & gdot for velocity vectors time td and ta

           f = 1 - xn^2*Cn/rd;
           f = 1 - yn/rd;
           g = tar - td - xn^3*Sn/sqrt(GM);
           g = A * sqrt(yn/GM);
           gdot = 1 - xn^2*Cn/ra;
           gdot =  1 - yn/ra;

           % Velocity vector at departure

           obj.ephemd(j).v = real((obj.ephema(j).r - f * obj.ephemd(j).r) / g) ;
           if obj.ni(j)>obj.nmax
               obj.ephemd(j).v = [Inf,Inf,Inf];
               obj.ephemd(j).V = Inf;
           else
               obj.ephemd(j).V = norm(obj.ephemd(j).v);           
           end


           % Velocity vector at arrival

           obj.ephema(j).v = real((- obj.ephemd(j).r + gdot * obj.ephema(j).r) / g) ;
           if obj.ni(j)>obj.nmax
               obj.ephema(j).v = [Inf,Inf,Inf];
               obj.ephema(j).V = Inf;
           else
               obj.ephema(j).V = norm(obj.ephema(j).v); 
           end
        end

       % Calculate Transfer Orbit

       obj.transfer_body(j).ephem0 = obj.ephemd(j);
       
       obj.transfer_body(j).ephemt = obj.ephemd(j);

       obj.transfer_body(j) = obj.transfer_body(j).calculate_orbit_from_ephem(td);
       obj.true_anom_dep(j) = obj.transfer_body(j).orbit.ta;

       obj.transfer_body(j).ephemt = obj.ephema(j);
       obj.transfer_body(j) = obj.transfer_body(j).calculate_orbit_from_ephem(tar);
       obj.true_anom_arr(j) = obj.transfer_body(j).orbit.ta;

    end
        
    if wayflag == 1
        if (dot([obj.ephemd(1).v(1),obj.ephemd(1).v(2)],[obj.bodyd.ephemt.v(1),obj.bodyd.ephemt.v(2)]) < 0 )
            if (dot([obj.ephemd(2).v(1),obj.ephemd(2).v(2)],[obj.bodyd.ephemt.v(1),obj.bodyd.ephemt.v(2)]) > 0 )
                obj.ephemd(1)=obj.ephemd(2);
                obj.ephema(1)=obj.ephema(2);
                obj.transfer_body(1)=obj.transfer_body(2);
                obj.true_anom_dep(1)=obj.true_anom_dep(2);
                obj.true_anom_arr(1)=obj.true_anom_arr(2);
            end
        elseif (dot([obj.ephemd(2).v(1),obj.ephemd(2).v(2)],[obj.bodyd.ephemt.v(1),obj.bodyd.ephemt.v(2)]) < 0 )
            obj.ephemd(2)=obj.ephemd(1);
            obj.ephema(2)=obj.ephema(1);
            obj.transfer_body(2)=obj.transfer_body(1);
            obj.true_anom_dep(2)=obj.true_anom_dep(1);
            obj.true_anom_arr(2)=obj.true_anom_arr(1);
        end
    elseif wayflag == -1
        if (dot([obj.ephemd(1).v(1),obj.ephemd(1).v(2)],[obj.bodyd.ephemt.v(1),obj.bodyd.ephemt.v(2)]) < 0 )
            if (dot([obj.ephemd(2).v(1),obj.ephemd(2).v(2)],[obj.bodyd.ephemt.v(1),obj.bodyd.ephemt.v(2)]) > 0 )
                obj.ephemd(2)=obj.ephemd(1);
                obj.ephema(2)=obj.ephema(1);
                obj.transfer_body(2)=obj.transfer_body(1);
                obj.true_anom_dep(2)=obj.true_anom_dep(1);
                obj.true_anom_arr(2)=obj.true_anom_arr(1);
            end
        elseif (dot([obj.ephemd(2).v(1),obj.ephemd(2).v(2)],[obj.bodyd.ephemt.v(1),obj.bodyd.ephemt.v(2)]) < 0 )
            obj.ephemd(1)=obj.ephemd(2);
            obj.ephema(1)=obj.ephema(2);
            obj.transfer_body(1)=obj.transfer_body(2);
            obj.true_anom_dep(1)=obj.true_anom_dep(2);
            obj.true_anom_arr(1)=obj.true_anom_arr(2);
        end        
    end
   return;
   


%%# FZERO solves the value of zn for which the yn is zero
function [ ZN, dzn ] = FZERO( z0, ztol, Maxit)
%# FZERO solves the value of zn for which the yn is zero
%#
%# INPUT:
%#
%# z0   : initial guess for auziliary variable zn
%# ztol : tolerance on zn for solution
%# A    : Universal Variable Formluation for A
%# Maxit: Maximum number of iterations
%#
%# OUTPUT :
%#
%# zn 	: value of zn for solution
%#      

	dzold = 0.0;
	ZN = z0;
    ZNOLD = ZN+ztol*10;
    count=0;

	reduce = 1.0;
	func1old=0.0;
	for i2=1:Maxit
		S1 = obj.bodya.Sz(ZN); 
		C1 = obj.bodya.Cz(ZN);

		func1 = rd + ra - A*(1.0 - ZN*S1) / sqrt(C1);

		grad = A / 4.0 * sqrt(C1) ;

		dz =  reduce*func1 / grad;
		ZN = ZN - dz;
        

		if ((func1old * func1) < 0.0&& count<5)
			reduce = reduce*0.1;
            ZN = ZNOLD;
            count=count+1;
        else
            count=0;
            reduce = 1.0;
        end
        if (abs(func1) < ztol)
            break;
        end
        dzold = dz;
        ZNOLD = ZN;
        func1old=func1;
    end
    dzn=abs(func1)/grad;
 %   i2
 %   func1

    return;
end


end %# Calculate_transfer

%%# deltaVd Calculate Departure DeltaVs
function obj = deltaVd(obj)
%# deltaVd Calculate Departure DeltaVs 
%#
%# INPUT:
%#
%# obj          : Current Transfer_orbit object in question
%#
%# OUTPUT:
%#
%# obj          : Transfer_orbit object with dVd's calculated for short way and long way
%#            
    for i=1:2
        if obj.ni(i)>=obj.nmax
            obj.dVd(i)= 1e50;
        else
            obj.dVd(i) = norm( obj.ephemd(i).v - obj.bodyd.ephemt.v );
        end
    end
    return;
end %# deltaVd

%%# deltaVa Calculates Arrival DeltaVs
function obj = deltaVa(obj)
%# deltaVa Calculate Arrival DeltaVs 
%#
%# INPUT:
%#
%# obj          : Current Transfer_orbit object in question
%#
%# OUTPUT:
%#
%# obj          : Transfer_orbit object with dVa's calculated for short way and long way
%#          

    for i=1:2
        if obj.ni(i)>=obj.nmax
            obj.dVa(i)= 1e50;
        else
            obj.dVa(i) = norm( obj.ephema(i).v - obj.bodya.ephemt.v );
        end
    end
    return;
end %# deltaVa

function obj = Calculate_Perihelion(obj)
    
    for j=1:2
        if (obj.true_anom_dep(j)*obj.true_anom_arr(j)<0)
            if (dot(obj.ephemd(j).r,obj.ephemd(j).v)<0)                
                obj.perihelion(j) = obj.transfer_body(j).orbit.a * ( 1 - obj.transfer_body(j).orbit.e );
            else
                obj.perihelion(j) = min(obj.ephemd(j).R,obj.ephema(j).R);
            end
        else
            if(obj.transfer_body(j).orbit.e<1)
                if(abs(obj.tar-obj.td)>obj.transfer_body(j).orbit.TP/2)
                    obj.perihelion(j) = obj.transfer_body(j).orbit.a * ( 1 - obj.transfer_body(j).orbit.e );
                else
                    obj.perihelion(j) = min(obj.ephemd(j).R,obj.ephema(j).R);
                end
            else
                obj.perihelion(j) = min(obj.ephemd(j).R,obj.ephema(j).R);
            end
        end
    end
    return;
end
                        
    
end %# methods

              
end

