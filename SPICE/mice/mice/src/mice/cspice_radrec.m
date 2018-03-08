%-Abstract
%
%   CSPICE_RADREC converts the right ascension, declination
%   coordinates of a location to rectangular (Cartesian)
%   coordinates.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY,
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING,
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      range   the value(s) describing the distance of the position
%              from the origin.
%
%              [1,n] = size(range); double = class(range)
%
%      ra      the value(s) describing the right ascension of the 
%              right ascension of the position:  the angular
%              distance measured toward the east from the prime meridian
%              to the meridian containing the input point. The direction
%              of increasing right ascension is from the +X axis towards
%              the +Y axis.
%
%              [1,n] = size(ra); double = class(ra)
%
%      dec     the value(s) describing the declination of the position as
%              measured in radians. This is the angular distance from the 
%              XY plane to the position.
%
%              [1,n] = size(dec); double = class(dec)
%
%              The range of `dec' is unrestricted.  Units are radians.
%
%   the call:
%
%      rectan = cspice_radrec( range, ra, dec)
%
%   returns:
%
%      rectan   the array(s) containing the rectangular coordinates of the
%               position(s).
%
%               [3,n] = size(rectan); double = class(rectan)
%
%               'rectan' returns with the same units associated with 'range'.
%
%               'rectan' returns with the same vectorization measure, N,
%                as 'range', 'ra', and 'dec'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a standard kernel set.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define a set of 15 right ascension-declination data sets
%      % pairs (in degrees) for the earth's pole and the array of
%      % corresponding ephemeris times J2000 TDB.
%      %
%      right_ascen = [ 180.003739,
%                      180.003205,
%                      180.002671,
%                      180.002137,
%                      180.001602,
%                      180.001068,
%                      180.000534,
%                      360.000000,
%                      359.999466,
%                      359.998932,
%                      359.998397,
%                      359.997863,
%                      359.997329,
%                      359.996795,
%                      359.996261 ];
%
%       dec        = [ 89.996751,
%                      89.997215,
%                      89.997679,
%                      89.998143,
%                      89.998608,
%                      89.999072,
%                      89.999536,
%                      90.000000,
%                      89.999536,
%                      89.999072,
%                      89.998607,
%                      89.998143,
%                      89.997679,
%                      89.997215,
%                      89.996751 ];
%
%       et         = [ -18408539.52023917,
%                      -15778739.49107254,
%                      -13148939.46190590,
%                      -10519139.43273926,
%                      -7889339.40357262,
%                      -5259539.37440598,
%                      -2629739.34523934,
%                       60.68392730,
%                       2629860.71309394,
%                       5259660.74226063,
%                       7889460.77142727,
%                       10519260.80059391,
%                       13149060.82976055,
%                       15778860.85892719,
%                       18408660.88809383 ];
%
%      %
%      % Create a 1xN array of radii, the length of a
%      % unit vector (1) the same size as the above arrays.
%      %
%      n_elements  = size(dec);
%      rad         = ones( 1,  n_elements(1) );
%      z_hat       = [0; 0; 1];
%
%      %
%      % Convert the RA/DEC values to radians.
%      %
%      right_ascen = right_ascen * cspice_rpd;
%      dec         = dec * cspice_rpd;
%
%      %
%      % Convert the angular description of the unit vectors to
%      % Cartesian.
%      %
%      pole        = cspice_radrec( rad, right_ascen', dec');
%
%      %
%      % Retrieve the transformation matrix from frames J2000 to
%      % IAU_EARTH.
%      %
%      mat         = cspice_pxform( 'J2000', 'IAU_EARTH', et');
%
%      %
%      % Rotate the 'pole' vector set into IAU_FRAME. All vectors
%      % should equal (to round-off) the z direction unit vector.
%      %
%
%      disp( ['      ET                x            y '   ...
%                            '          z      Angular diff'] )
%      disp( [' ________________  __________  __________' ...
%                              '  __________ ______________'] )
%
%      for i =1:15
%         z = mat(:,:,i) * pole(:,i);
%
%         %
%         % Output the ephemeris time, the pole vector, and the angular
%         % separation between the calculated and the expected pole vectors.
%         %
%         txt = sprintf( '%18.8f %11.8f %11.8f %11.8f %11.8e', ...
%                        et(i), z, cspice_vsep(z,z_hat) );
%         disp(txt)
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%            ET                x            y           z      Angular diff
%       ________________  __________  __________  __________ ______________
%
%      -18408539.52023917  0.00000000 -0.00000000  1.00000000 2.72257100e-09
%      -15778739.49107254  0.00000000 -0.00000000  1.00000000 1.86400956e-10
%      -13148939.46190590 -0.00000000  0.00000000  1.00000000 3.09537269e-09
%      -10519139.43273926  0.00000000 -0.00000001  1.00000000 6.00434486e-09
%       -7889339.40357262  0.00000000 -0.00000001  1.00000000 8.53997578e-09
%       -5259539.37440598 -0.00000000  0.00000001  1.00000000 5.63100382e-09
%       -2629739.34523934 -0.00000000 -0.00000000  1.00000000 2.72203209e-09
%             60.68392730 -0.00000000 -0.00000000  1.00000000 1.86939958e-10
%        2629860.71309394  0.00000000  0.00000000  1.00000000 3.09591191e-09
%        5259660.74226063 -0.00000000 -0.00000001  1.00000000 6.00488364e-09
%        7889460.77142727 -0.00000000 -0.00000001  1.00000000 8.53943655e-09
%       10519260.80059391  0.00000000  0.00000000  1.00000000 5.63046483e-09
%       13149060.82976055 -0.00000000 -0.00000000  1.00000000 2.72149287e-09
%       15778860.85892719 -0.00000000 -0.00000000  1.00000000 1.87478860e-10
%       18408660.88809383  0.00000000  0.00000000  1.00000000 3.09645104e-09
%
%   The angular deviation between the calculated pole vector and the expected
%   measures as ~10**-9.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine radrec_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.2, 07-JAN-2015, EDW (JPL)
%
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   range ra and dec to rectangular coordinates
%   right_ascension and declination to rectangular
%
%-&

function [rectan] = cspice_radrec( range, ra, dec)

   switch nargin
      case 3

         range = zzmice_dp(range);
         ra    = zzmice_dp(ra);
         dec   = zzmice_dp(dec);

      otherwise
         error ( 'Usage: [_rectan(3)_] = cspice_radrec( _range_, _ra_, _dec_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice('radrec_c', range, ra, dec);
   catch
      rethrow(lasterror)
   end

