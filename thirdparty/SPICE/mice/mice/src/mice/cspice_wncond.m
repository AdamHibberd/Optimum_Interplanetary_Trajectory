%-Abstract
%
%   CSPICE_WNCOND contracts each of the intervals of a double
%   precision window.
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
%      left     amount to add to each left interval endpoint in the input
%               'window'.
%
%               [1,1] = size(left); double = class(left)
%
%      right    amount to subtract to each right interval endpoint in the input
%               'window'.
%
%               [1,1] = size(right); double = class(right)
%
%      window   SPICE window containing zero or more intervals.
%
%               [2m,1] = size(window); double = class(window)
%
%   the call:
%
%      window_f = cspice_wncond( left, right, window)
%
%   returns:
%
%      window_f   SPICE window containing zero or more intervals, representing
%                 the original 'window' with each of its intervals contracted
%                 by 'left' units on the left and 'right' units on the right
%
%                 [2n,1] = size(window_f); double = class(window_f)
%
%                 'window_f' can overwrite 'window'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Let 'window' contain the intervals
%      %
%      window = [ [ 1; 3 ]; [ 7; 11 ]; [ 23; 27 ]; [ 29; 29 ]; ];
%
%      %
%      % Apply the following series of calls:
%      %
%      window = cspice_wncond(  2,  1, window )
%      window = cspice_wncond( -2,  2, window )
%      window = cspice_wncond( -2, -1, window )
%
%   MATLAB outputs:
%
%      window =
%
%           9
%          10
%          25
%          26
%
%      Representing the intervals:
%
%         [ 9, 10 ]  [ 25, 26 ]
%
%      window =
%
%           7
%           8
%          23
%          24
%
%      Representing the intervals:
%
%         [ 7,  8 ]  [ 23, 24 ]
%
%      window =
%
%           5
%           9
%          21
%          25
%
%      Representing the intervals:
%
%         [ 5,  9 ]  [ 21, 25 ]
%
%
%   Note that intervals may be "contracted" by negative amounts.
%   In the example above, the second call shifts each interval to
%   the left, while the third call undoes the effect of the first
%   call (without restoring the destroyed intervals).
%
%   Note also that the third call is exactly equivalent to the
%   call:
%
%      cspice_wnexpd( 2, 1, window )
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wncond_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 24-JUL-2007, EDW (JPL)
%
%-Index_Entries
%
%   contract the intervals of a d.p. window
%
%-&

function [window_f] = cspice_wncond( left, right, window)

   switch nargin

      case 3

         left   = zzmice_dp(left);
         right  = zzmice_dp(right);
         window = zzmice_win(window);

      otherwise

         error ( 'Usage: [window_f] = cspice_wncond( left, right, window)' )

   end

   %
   % Please note, this call does not require addition of space for the window
   % control segment as needed by other windows interfaces. The interface
   % copies the data in 'window' to a work variable rather than directly
   % pass 'window' into a CSPICE call.
   %
   % Call the wnexpd_c interface with negated values of 'left' and 'right'.
   %
   try
      [window_f] = mice('wnexpd_c', -left, -right, window );
   catch
      rethrow(lasterror)
   end




