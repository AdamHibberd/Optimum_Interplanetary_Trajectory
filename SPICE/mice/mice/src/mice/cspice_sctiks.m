%-Abstract
%
%   CSPICE_SCTIKS convert a spacecraft clock format string to
%   number of 'ticks'.
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
%      sc       the NAIF ID of the spacecraft clock whose time is 
%               to encode.
%
%               [1,1] = size(sc); int32 = class(sc)
%
%      clkstr   the scalar string or N-vector representation of the
%               'sc' spacecraft's clock time, WITHOUT PARTITION NUMBER.
%
%               [n,c1] = size(clkstr); char = class(clkstr)
%
%                  or
%
%               [1,n] = size(clkstr); cell = class(clkstr)
%
%   the call:
%
%      ticks = cspice_sctiks( sc, clkstr )
%
%   returns:
%
%      ticks   the tick values(s) represented by the spacecraft clock 
%              string 'clkstr'.
%
%              [1,n] = size(ticks); double = class(ticks)
%
%              'ticks' returns with the same vectorization measure, N,
%              as 'clkstr'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%
%   MATLAB outputs:
%
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine sctiks_c.
%
%   MICE.REQ
%   SCLK.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 06-JAN-2015, EDW (JPL)
%
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 07-JUN-2006, EDW (JPL)
%
%-Index_Entries
%
%   convert spacecraft_clock string to ticks
%
%-&

function [ticks] = cspice_sctiks(sc, clkstr)

   switch nargin
      case 2

         sc     = zzmice_int(sc);
         clkstr = zzmice_str(clkstr);

      otherwise

         error ( 'Usage: [_ticks_] = cspice_sctiks(sc, _`clkstr`_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [ticks] = mice('sctiks_c', sc, clkstr);
   catch
      rethrow(lasterror)
   end





