%-Abstract
%
%   MICE_SUBPNT computes the rectangular coordinates of the
%   sub-observer point on a target body at a specified epoch,
%   optionally corrected for light time and stellar aberration.
%
%   This routine supersedes mice_subpt, which does not have an input
%   argument for the target body-fixed frame name.
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
%      method   a string providing parameters defining
%               the computation method to use.
%
%               [1,c1] = size(method); char = class(method)
%
%                  or
%
%               [1,1] = size(method); cell = class(method)
%
%               The supported values of 'method' are listed below.
%               Please note that the colon is a required delimiter;
%               using a blank will not work.
%
%                     'Near point: ellipsoid'   The sub-observer point
%                                               computation uses a
%                                               triaxial ellipsoid to
%                                               model the surface of the
%                                               target body. The
%                                               sub-observer point is
%                                               defined as the nearest
%                                               point on the target
%                                               relative to the observer.
%
%                     'Intercept: ellipsoid'    The sub-observer point
%                                               computation uses a
%                                               triaxial ellipsoid to
%                                               model the surface of the
%                                               target body. The
%                                               sub-observer point is
%                                               defined as the target
%                                               surface intercept of the
%                                               line containing the
%                                               observer and the
%                                               target's center.
%
%               Neither case nor white space are significant in
%               'method'. For example, the string
%
%                    ' nearpoint:ELLIPSOID '
%
%               is valid.
%
%      target   the name of the target body. The target
%               body is an ephemeris object (its trajectory is given by
%               SPK data), and is an extended object.
%
%               [1,c2] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               The string 'target' is case-insensitive, and leading
%               and trailing blanks in 'target' are not significant.
%               Optionally, you may supply a string containing the
%               integer ID code for the object. For example both
%               'MOON' and '301' are legitimate strings that indicate
%               the moon is the target body.
%
%               When the target body's surface is represented by a
%               tri-axial ellipsoid, this routine assumes that a
%               kernel variable representing the ellipsoid's radii is
%               present in the kernel pool. Normally the kernel
%               variable would be defined by loading a PCK file.
%
%
%      et       the epoch(s), expressed as seconds past
%               J2000 TDB, of the observer: 'et' is
%               the epoch at which the observer's state is computed.
%
%               [1,n] = size(et); double = class(et)
%
%               When aberration corrections are not used, 'et' is also
%               the epoch at which the position and orientation of
%               the target body are computed.
%
%               When aberration corrections are used, 'et' is the epoch
%               at which the observer's state relative to the solar
%               system barycenter is computed; in this case the
%               position and orientation of the target body are
%               computed at et-lt or et+lt, where 'lt' is the one-way
%               light time between the sub-observer point and the
%               observer, and the sign applied to 'lt' depends on the
%               selected correction. See the description of 'abcorr'
%               below for details.
%
%      fixref   the name of the body-fixed, body-centered
%               reference frame associated with the target body.
%               The output sub-observer point 'spoint' will be
%               expressed relative to this reference frame.
%
%               [1,c3] = size(fixref); char = class(fixref)
%
%                  or
%
%               [1,1] = size(fixref); cell = class(fixref)
%
%      abcorr   the aberration correction to apply
%               when computing the observer-target state and the
%               orientation of the target body.
%
%               [1,c4] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               For remote sensing applications, where the apparent
%               sub-observer point seen by the observer is desired,
%               normally either of the corrections
%
%                     'LT+S'
%                     'CN+S'
%
%               should be used. These and the other supported options
%               are described below. 'abcorr' may be any of the
%               following:
%
%                     'NONE'     Apply no correction. Return the
%                                geometric sub-observer point on the
%                                target body.
%
%               Let 'lt' represent the one-way light time between the
%               observer and the sub-observer point (note: NOT
%               between the observer and the target body's center).
%               The following values of 'abcorr' apply to the
%               "reception" case in which photons depart from the
%               sub-observer point's location at the light-time
%               corrected epoch et-lt and *arrive* at the observer's
%               location at 'et':
%
%                     'LT'       Correct for one-way light time (also
%                                called "planetary aberration") using a
%                                Newtonian formulation. This correction
%                                yields the location of sub-observer
%                                point at the moment it emitted photons
%                                arriving at the observer at 'et'.
%
%                                The light time correction uses an
%                                iterative solution of the light time
%                                equation. The solution invoked by the
%                                'LT' option uses one iteration.
%
%                                Both the target position as seen by the
%                                observer, and rotation of the target
%                                body, are corrected for light time.
%
%                     'LT+S'     Correct for one-way light time and
%                                stellar aberration using a Newtonian
%                                formulation. This option modifies the
%                                state obtained with the 'LT' option to
%                                account for the observer's velocity
%                                relative to the solar system
%                                barycenter. The result is the apparent
%                                sub-observer point as seen by the
%                                observer.
%
%                     'CN'       Converged Newtonian light time
%                                correction. In solving the light time
%                                equation, the 'CN' correction iterates
%                                until the solution converges. Both the
%                                position and rotation of the target
%                                body are corrected for light time.
%
%                     'CN+S'     Converged Newtonian light time and
%                                stellar aberration corrections. This
%                                option produces a solution that is at
%                                least as accurate at that obtainable
%                                with the 'LT+S' option. Whether the 'CN+S'
%                                solution is substantially more accurate
%                                depends on the geometry of the
%                                participating objects and on the
%                                accuracy of the input data. In all
%                                cases this routine will execute more
%                                slowly when a converged solution is
%                                computed.
%
%      obsrvr   the scalar string name of the observing body. The
%               observing body is an ephemeris object: it typically
%               is a spacecraft, the earth, or a surface point on the
%               earth. 'obsrvr' is case-insensitive, and leading and
%               'obsrvr' are not significant. Optionally, you may
%               trailing blanks in supply a string containing the integer
%               ID code for the object. For example both 'MOON' and '301'
%               are legitimate strings that indicate the Moon is the
%               observer.
%
%               [1,c5] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%   the call:
%
%      [subpnt] = mice_subpnt( method, target, et, fixref, ...
%                              abcorr, obsrvr)
%
%   returns:
%
%      subpnt  the structure(s) containing the results of the calculation.
%
%              [1,n] = size(subpnt); struct = class(subpnt)
%
%              Each structure consists of the fields:
%
%              'spoint'   the array defining the sub-observer point on the
%                         target body.
%
%                         [3,1] = size(subpnt(i).spoint)
%                         double = class(subpnt(i).spoint)
%
%                         The sub-observer point is defined either as the point
%                         on the target body that is closest to the observer,
%                         or the target surface intercept of the line from the
%                         observer to the target's center; the input argument
%                         'method' selects the definition to be used.
%
%                         'spoint' is expressed in Cartesian coordinates,
%                         relative to the body-fixed target frame designated by
%                         'fixref'. The body-fixed target frame is evaluated at
%                         the sub-observer epoch 'trgepc'
%                         (see description below).
%
%                         When light time correction is used, the duration of
%                         light travel between 'spoint' to the observer is
%                         considered to be the one way light time.
%
%                         When aberration corrections are used, 'spoint'
%                         is computed using target body position and
%                         orientation that have been adjusted for the
%                         corrections applicable to 'spoint' itself
%                         rather than to the target body's center. In
%                         particular, if the stellar aberration correction
%                         applicable to 'spoint' is represented by a shift
%                         vector 's', then the light-time corrected
%                         position of the target is shifted by 's' before
%                         the sub-observer point is computed.
%
%                         The components of 'spoint' have units of km.
%
%              'trgepc'   the "sub-observer point epoch." 'trgepc' is defined
%                         as follows: letting 'lt' be the one-way 'trgepc'
%                         is the epoch et-lt, et+lt, or 'et' depending on
%                         whether the requested aberration correction is,
%                         respectively, for received radiation, transmitted
%                         radiation, or omitted. 'lt' is computed using
%                         the method indicated by 'abcorr'.
%
%                         [1,1] = size(subpnt(i).trgepc)
%                         double = class(subpnt(i).trgepc)
%
%                         'trgepc' is expressed as seconds past J2000 TDB.
%
%              'srfvec'   the array defining the position vector from
%                         the observer at 'et' to 'spoint'. 'srfvec' is
%                         expressed in the target body-fixed  reference frame
%                         designated by 'fixref', evaluated at  'trgepc'.
%
%                         [3,1] = size(subpnt(i).srfvec)
%                         double = class(subpnt(i).srfvec)
%
%                         The components of 'srfvec' are given in units of km.
%
%                         One can use the CSPICE function vnorm_c to obtain the
%                         distance between the observer and 'spoint':
%
%                            dist = norm( subpnt(i).srfvec )
%
%                         The observer's position 'obspos', relative to the
%                         target body's center, where the center's position is
%                         corrected for aberration effects as indicated by
%                         'abcorr', can be computed with:
%
%                            obspos =  subpnt(i).spoint - subpnt(i).srfvec
%
%                         To transform the vector 'srfvec' to a time-dependent
%                         reference frame 'ref' at 'et', a sequence of two
%                         frame transformations is required. For example, let
%                         'mfix' and 'mref' be 3x3 matrices respectively
%                         describing the target body-fixed to inertial frame
%                         transformation at 'trgepc' and the inertial to
%                         (time-dependent frame) 'ref' transformation at 'et',
%                         and let 'xform' be the 3x3 matrix representing the
%                         composition of 'mref' with 'mfix'. Then 'srfvec' can
%                         be transformed to the result 'refvec' as follows:
%
%                            mfix   = cspice_pxform( fixref,       ...
%                                                    'j2000',      ...
%                                                    subpnt(i).trgepc )
%                            mref   = cspice_pxform( 'j2000', ref, et )
%                            xform  = mref * mfix
%                            refvec = xform * subpnt(i).srfvec
%
%      'subpnt' return with the same vectorization measure, N, as 'et'.
%
%      Note, if needed, the user can extract the field data from vectorized
%      'spoint' structures into separate arrays.
%
%      Extract the 'spoint' field data to a 3X1 array 'spoint':
%
%         spoint = reshape( [subpnt(:).spoint], 3, [] )
%
%      Extract the 'trgepc' field data to a scalar 'trgepc':
%
%         trgepc = reshape( [subpnt(:).trgepc], 1, [] )
%
%      Extract the 'spoint' field data to a 3X1 array 'spoint':
%
%         spoint = reshape( [subpnt(:).spoint], 3, [] )
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load kernel files via the meta-kernel.
%      %
%      cspice_furnsh( '/kernels/standard.tm' );
%
%      %
%      % Convert the UTC request time to ET (seconds past
%      % J2000, TDB).
%      %
%      et0 = cspice_str2et( '2008 aug 11 00:00:00' );
%
%      %
%      % Create a vector of times. The code will also run for 'et'
%      % a scalar.
%      %
%      et = [0:10]*cspice_spd + et0;
%
%      %
%      % Look up the target body's radii. We'll use these to
%      % convert Cartesian to planetodetic coordinates. Use
%      % the radii to compute the flattening coefficient of
%      % the reference ellipsoid.
%      %
%      radii  = cspice_bodvrd( 'MARS', 'RADII', 3 );
%
%      %
%      % Let RE and RP be, respectively, the equatorial and
%      % polar radii of the target.
%      %
%      re = radii(1);
%      rp = radii(3);
%      f = ( re-rp)/re;
%
%      %
%      % Compute sub-observer point using light time and stellar
%      % aberration corrections. Use the "target surface intercept"
%      % definition of the sub-observer point on the first loop
%      % iteration, and use the "near point" definition on the
%      % second.
%      %
%
%      method = { 'Intercept:  ellipsoid', 'Near point: ellipsoid' };
%
%      for i=1:2
%
%         subpnt = mice_subpnt( method(i), 'MARS', et, ...
%                              'IAU_MARS', 'LT+S', 'EARTH' );
%
%         N = length(subpnt);
%
%         %
%         % Expand the embedded data arrays to properly shaped
%         % generic arrays.
%         %
%         spoint   = reshape( [subpnt.spoint], 3, [] );
%         trgepc   = reshape( [subpnt.trgepc], 1, [] );
%         srfvec   = reshape( [subpnt.srfvec], 3, [] );
%
%         %
%         % Convert the sub-observer point's rectangular coordinates
%         % to planetographic longitude, latitude and altitude.
%         % Convert radians to degrees.
%         %
%         [ spglon, spglat, spgalt] = cspice_recpgr( 'mars', spoint, ...
%                                                     re,    f );
%
%         spglon = spglon * cspice_dpr;
%         spglat = spglat * cspice_dpr;
%
%         %
%         % Convert sub-observer point's rectangular coordinates to
%         % planetocentric radius, longitude, and latitude. Convert
%         % radians to degrees.
%         %
%         [ spcrad, spclon, spclat ] =cspice_reclat( spoint ) ;
%
%         spclon = spclon * cspice_dpr;
%         spclat = spclat * cspice_dpr;
%
%         %
%         % Compute the observer's position relative to the center of the
%         % target, where the center's location has been adjusted using
%         % the aberration corrections applicable to the sub-point.
%         % Express the observer's location in geodetic coordinates.
%         %
%         obspos = spoint - srfvec;
%
%         [ opglon, opglat, opgalt] = cspice_recpgr( 'mars', obspos, re, f );
%
%         opglon = opglon * cspice_dpr;
%         opglat = opglat * cspice_dpr;
%
%         %
%         % Convert the observer's rectangular coordinates to planetodetic
%         % longitude, latitude and altitude. Convert radians to degrees.
%         %
%         [opcrad, opclon, opclat] = cspice_reclat( obspos ) ;
%
%         opclon = opclon * cspice_dpr;
%         opclat = opclat * cspice_dpr;
%
%         utcstr = cspice_et2utc( et, 'C', 6);
%
%         for j=1:N
%
%           fprintf( 'Computational Method %s\n\n', char(method(i)) )
%
%           fprintf( 'Time (UTC):                          %s\n',  ...
%                                                        utcstr(j,:) )
%           fprintf(                                                  ...
%           'Observer altitude                      (km) = %21.9f\n', ...
%                                                        opgalt(j) )
%
%           fprintf(                                                  ...
%           'Length of SRFVEC                       (km) = %21.9f\n', ...
%                                               norm(srfvec(:,j))  )
%
%           fprintf(                                                  ...
%           'Sub-observer point altitude            (km) = %21.9f\n', ...
%                                                        spgalt(j) )
%
%           fprintf(                                                  ...
%           'Sub-observer planetographic longitude (deg) = %21.9f\n', ...
%                                                        spglon(j) )
%
%           fprintf(                                                  ...
%           'Observer planetographic longitude     (deg) = %21.9f\n', ...
%                                                        opglon(j) )
%
%           fprintf(                                                  ...
%           'Sub-observer planetographic latitude  (deg) = %21.9f\n', ...
%                                                        spglat(j) )
%
%           fprintf(                                                  ...
%           'Observer planetographic latitude      (deg) = %21.9f\n', ...
%                                                        opglat(j) )
%
%           fprintf(                                                  ...
%           'Sub-observer planetocentric longitude (deg) = %21.9f\n', ...
%                                                        spclon(j) )
%
%           fprintf(                                                  ...
%           'Observer planetocentric longitude     (deg) = %21.9f\n', ...
%                                                        opclon(j) )
%
%           fprintf(                                                  ...
%           'Sub-observer planetocentric latitude  (deg) = %21.9f\n', ...
%                                                        spclat(j) )
%
%           fprintf(                                                  ...
%           'Observer planetocentric latitude      (deg) = %21.9f\n', ...
%                                                        opclat(j) )
%
%           fprintf( '\n')
%
%         end
%
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
%      Computational Method Intercept:  ellipsoid
%
%      Time (UTC):                          2008 AUG 11 00:00:00.000000
%      Observer altitude                      (km) =   349199089.540775955
%      Length of SRFVEC                       (km) =   349199089.577471614
%      Sub-observer point altitude            (km) =           0.000000000
%      Sub-observer planetographic longitude (deg) =         199.302305032
%      Observer planetographic longitude     (deg) =         199.302305032
%      Sub-observer planetographic latitude  (deg) =          26.262401237
%      Observer planetographic latitude      (deg) =          25.994936751
%      Sub-observer planetocentric longitude (deg) =         160.697694968
%      Observer planetocentric longitude     (deg) =         160.697694968
%      Sub-observer planetocentric latitude  (deg) =          25.994934171
%      Observer planetocentric latitude      (deg) =          25.994934171
%
%         ...
%
%      Computational Method Intercept:  ellipsoid
%
%      Time (UTC):                          2008 AUG 21 00:00:00.000212
%      Observer altitude                      (km) =   355208039.374441564
%      Length of SRFVEC                       (km) =   355208039.410808742
%      Sub-observer point altitude            (km) =           0.000000000
%      Sub-observer planetographic longitude (deg) =         101.105921254
%      Observer planetographic longitude     (deg) =         101.105921254
%      Sub-observer planetographic latitude  (deg) =          26.097147650
%      Observer planetographic latitude      (deg) =          25.830877512
%      Sub-observer planetocentric longitude (deg) =        -101.105921254
%      Observer planetocentric longitude     (deg) =        -101.105921254
%      Sub-observer planetocentric latitude  (deg) =          25.830874987
%      Observer planetocentric latitude      (deg) =          25.830874987
%
%      Computational Method Near point: ellipsoid
%
%      Time (UTC):                          2008 AUG 11 00:00:00.000000
%      Observer altitude                      (km) =   349199089.540767729
%      Length of SRFVEC                       (km) =   349199089.540767789
%      Sub-observer point altitude            (km) =          -0.000000000
%      Sub-observer planetographic longitude (deg) =         199.302305032
%      Observer planetographic longitude     (deg) =         199.302305032
%      Sub-observer planetographic latitude  (deg) =          25.994936751
%      Observer planetographic latitude      (deg) =          25.994936751
%      Sub-observer planetocentric longitude (deg) =         160.697694968
%      Observer planetocentric longitude     (deg) =         160.697694968
%      Sub-observer planetocentric latitude  (deg) =          25.729407227
%      Observer planetocentric latitude      (deg) =          25.994934171
%
%         ...
%
%      Computational Method Near point: ellipsoid
%
%      Time (UTC):                          2008 AUG 21 00:00:00.000212
%      Observer altitude                      (km) =   355208039.374358773
%      Length of SRFVEC                       (km) =   355208039.374358773
%      Sub-observer point altitude            (km) =          -0.000000000
%      Sub-observer planetographic longitude (deg) =         101.105921255
%      Observer planetographic longitude     (deg) =         101.105921255
%      Sub-observer planetographic latitude  (deg) =          25.830877512
%      Observer planetographic latitude      (deg) =          25.830877512
%      Sub-observer planetocentric longitude (deg) =        -101.105921255
%      Observer planetocentric longitude     (deg) =        -101.105921255
%      Sub-observer planetocentric latitude  (deg) =          25.566547790
%      Observer planetocentric latitude      (deg) =          25.830874987
%
%-Particulars
%
%   A sister version of this routine exists named cspice_subpnt that returns
%   the structure field data as separate arguments.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine subpnt_c.
%
%   MICE.REQ
%   FRAMES.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.1.0, 12-JAN-2015, EDW (JPL)
%
%       Vectorized interface on input 'et'.
%
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.1, 12-MAY-2009, EDW (JPL)
%
%       Corrected type in I/O call description. The call description
%       lacked the 'fixref' argument.
%
%   -Mice Version 1.0.0, 30-JAN-2008, EDW (JPL)
%
%-Index_Entries
%
%   find sub-observer point on target body
%   find nearest point to observer on target body
%
%-&

function [subpnt] = mice_subpnt( method, target, et, fixref, abcorr, obsrvr )

   switch nargin
      case 6

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);

      otherwise
         error ( ['Usage: [_subnt_] = '                          ...
                  'mice_subpnt( `method`, `target`, _et_,'       ...
                  ' `fixref`, `abcorr`, `obsrvr`)'])
   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [subpnt] = mice('subpnt_s', method, target, et, fixref, abcorr, obsrvr);
   catch
      rethrow(lasterror)
   end


