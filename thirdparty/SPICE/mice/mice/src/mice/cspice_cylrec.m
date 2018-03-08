%-Abstract
%
%   CSPICE_CYLREC converts cylindrical coordinates to rectangular
%   (Cartesian) coordinates.
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
%      r      the value(s) describing the distance of the point of
%             interest from z axis.
%
%             [1,n] = size(r); double = class(r)
%
%      lonc   the value(s) describing the cylindrical angle of the point of
%             interest from the XZ plane measured in radians.
%
%             [1,n] = size(lonc); double = class(lonc)
%
%      z      the value(s) describing the height of the point above
%             the XY plane.
%
%             [1,n] = size(z); double = class(z)
%
%   the call:
%
%      rectan = cspice_cylrec( r, lonc, z)
%
%   returns:
%
%      rectan   the array(s) containing the rectangular coordinates of the
%               position or set of positions
%
%               [3,n] = size(rectan); double = class(rectan)
%
%               The argument 'rectan' returns in the same units associated
%               with 'r' and 'z'.
%
%               'rectan' returns with the same vectorization measure (N) as
%               'r', 'lonc', and 'z'
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example (1):
%
%      %
%      % Load an SPK, leapseconds, and PCK kernel set.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Create a vector of scalar times.
%      %
%      et = [0:2]*2.*cspice_spd;
%
%      %
%      % Retrieve the position of the moon seen from earth at 'et'
%      % in the J2000 frame without aberration correction.
%      %
%      [pos, et] = cspice_spkpos( 'MOON', et, 'J2000', 'NONE', 'EARTH' );
%
%      %
%      % Convert the array of position vectors 'pos' to cylindrical
%      % coordinates.
%      %
%      [r, lonc, z] = cspice_reccyl(pos);
%
%      %
%      % Convert the cylindrical to rectangular.
%      %
%      [rectan] = cspice_cylrec(r, lonc, z);
%
%      %
%      % Calculate the relative error against the original position
%      % vectors.
%      %
%      (rectan-pos) ./ pos
%
%      end
%
%   MATLAB outputs:
%
%      1.0e-13 *
%
%      0.00199609007208                  0  -0.25513381329527
%     -0.00218237675815                  0  -0.00153127196389
%                     0                  0                  0
%
%       The relative error between the original array of position vectors
%       and those resulting from the coordinate conversions
%       has magnitude on the order of 10^(-13).  A numerical
%       demonstration of equality.
%
%   Example (2):
%
%      %
%      % Define six sets of cylindrical coordinates, 'lonc' expressed
%      % in degrees - converted to radians by use of cspice_rpd.
%      %
%      r     = [ 1.,  1.,   1.,   1.,   0.,  0. ];
%      lonc  = [ 0., 90., 180., 180., 180., 33. ] * cspice_rpd;
%      z     = [ 0.,  0.,   1.,  -1.,   1.,  0. ];
%
%      %
%      % ...convert the cylindrical coordinates to rectangular coordinates
%      %
%      [rec] = cspice_cylrec(r, lonc, z);
%
%      %
%      % ...convert angular measure to degrees.
%      %
%      lonc = lonc * cspice_dpr;
%
%      disp('     r         lonc        z           x         y           z   ')
%      disp('  --------   --------   --------   --------   --------   --------')
%
%      %
%      % Create an array of values for output.
%      %
%      output = [ r; lonc; z; rec(1,:); rec(2,:); rec(3,:) ];
%
%      txt = sprintf( '%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', output );
%      disp( txt )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%        r         lonc        z           x         y           z
%     --------   --------   --------   --------   --------   --------
%       1.0000     0.0000     0.0000     1.0000     0.0000     0.0000
%       1.0000    90.0000     0.0000     0.0000     1.0000     0.0000
%       1.0000   180.0000     1.0000    -1.0000     0.0000     1.0000
%       1.0000   180.0000    -1.0000    -1.0000     0.0000    -1.0000
%       0.0000   180.0000     1.0000    -0.0000     0.0000     1.0000
%       0.0000    33.0000     0.0000     0.0000     0.0000     0.0000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine cylrec_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 30-OCT-2014, EDW (JPL)
%
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   cylindrical to rectangular
%
%-&

function [rectan] = cspice_cylrec(r, lonc, z)

   switch nargin
      case 3

         r    = zzmice_dp(r);
         lonc = zzmice_dp(lonc);
         z    = zzmice_dp(z);

      otherwise

         error ( 'Usage: [_rectan(3)_] = cspice_cylrec(_r_, _lonc_, _z_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice('cylrec_c', r, lonc, z);
   catch
      rethrow(lasterror)
   end

