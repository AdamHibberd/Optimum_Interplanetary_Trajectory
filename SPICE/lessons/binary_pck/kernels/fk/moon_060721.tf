
KPL/FK


   SPICE Lunar Reference Frame Specification Kernel
   =====================================================================

   Original file name:                   moon_060721.tf
   Creation date:                        2006 July 21 16:52
   Revision date:                        2007 February 12 03:14
   Created by:                           Nat Bachman  (NAIF/JPL)

   Revision description:

      The comments in this kernel have been updated to refer to the
      current lunar binary PCK.  Also, this kernel now refers to 
      the lunar frame association kernels which direct the SPICE 
      system to associate with the Moon either the MOON_ME or MOON_PA 
      reference frames.

      The data (that is, the "keyword = value" assignments) in this
      kernel are unchanged.
 

   Frames Defined by this Kernel
   =====================================================================

   Frame Name       Relative to        Type   Frame ID
   --------------   -----------------  -----  --------
   MOON_PA          MOON_PA_DE403      FIXED  31000
   MOON_ME          MOON_ME_DE403      FIXED  31001
   MOON_PA_DE403    ICRF/J2000         PCK    31002
   MOON_ME_DE403    MOON_PA_DE403      FIXED  31003
        

   Introduction
   =====================================================================

   This kernel specifies lunar body-fixed reference frames for use by
   SPICE-based application software.  These reference frames are
   associated with high-accuracy lunar orientation data provided by the
   JPL Solar System Dynamics Group's planetary ephemerides (both
   trajectory and lunar orientation data are stored in these ephemeris
   files).  These ephemerides have names of the form DE-nnn (DE stands
   for "developmental ephemeris").

   The frames specified by this kernel are realizations of two different
   lunar reference systems:

      Principal axes (PA) system 
      --------------------------
      The axes of this system are defined by the principal axes of the
      Moon.  Due to the nature of the Moon's orbit and
      rotation, the Z axis of this system does not coincide with the
      Moon's mean spin axis, nor does the X axis coincide with the mean
      direction to the center of the Earth (in contrast with the ME
      system defined below).
 
      Lunar principal axes frames realizing the lunar PA system and
      specified by this kernel are associated with JPL planetary
      ephemerides.  Each new JPL planetary ephemeris can (but does not
      necessarily) define a new realization of the lunar principal axes
      system.  Coordinates of lunar surface features expressed in lunar
      PA frames can change slightly from one lunar ephemeris version to
      the next.
 

      Mean Earth/polar axis (ME) system
      ---------------------------------
      The Lunar mean Earth/axis system is a lunar body-fixed reference
      system used in the IAU/IAG Working Group Report [2] to describe
      the orientation of the Moon relative to the ICRF frame. The +Z
      axis of this system is aligned with the mean lunar north pole,
      while the prime meridian contains the the mean Earth direction.
 
      The mean directions used to define the axes of a mean Earth/polar
      axis reference frame realizing the lunar ME system and specified
      by this kernel are associated with a given JPL planetary
      ephemeris version.  The rotation between the mean Earth frame for
      a given ephemeris version and the associated principal axes frame
      is given by a constant matrix (see [1]).


   For each JPL planetary ephemeris (DE), this kernel includes
   specifications of the corresponding principal axes and mean Earth/
   polar axis frames.  The names of these frames have the form

      MOON_PA_DEnnn
     
   and

      MOON_ME_DEnnn

   respectively, where nnn is the version number of the DE. This set of
   DE-dependent frame specifications will grow over time; frame
   specifications pertaining to older DEs will be retained in later
   versions of this frame kernel.

   For each of the two reference systems, there is a corresponding
   "generic" frame specification:  these generic frames are simply
   aliases for the PA and ME frames associated with the latest DE.  The
   generic frame names are

      MOON_PA
      MOON_ME

   These generic frame names are provided to enable SPICE-based
   applications to refer to the latest DE-based (or other) lunar
   rotation data without requiring code modifications as new data
   become available. SPICE users may, if they wish, modify this kernel
   to assign these frame aliases to other frames than those selected
   here, for example, older DE-based frames.  NAIF recommends that, if
   this frame kernel is modified, the name of this file also be changed to
   avoid confusion.
 

   Comparison of PA and ME frames
   ------------------------------

   The rotation between the mean Earth frame for a given DE and the
   associated principal axes frame for the same DE is given by a constant
   matrix (see [1]).  For DE-403, the rotation angle of this matrix is
   approximately 0.028241 degrees; this is equivalent to approximately 860 m
   when expressed as a displacement along a great circle on the Moon's
   surface.


   Comparison of DE-based and IAU/IAG report rotation data
   -------------------------------------------------------

   Within the SPICE system, the lunar ME frame specified by the
   rotational elements from the IAU/IAG Working Group report [2] is
   given the name IAU_MOON; the data defining this frame are provided
   in a generic text PCK.

   The orientation of the lunar ME frame obtained by applying the
   DE-based PA-to-ME rotation described above to the DE-based lunar
   libration data does not agree closely with the lunar ME frame
   orientation given by the rotational elements from the IAU/IAG
   Working Group report (that is, the IAU_MOON frame). The difference
   is due to truncation of the libration series used in the report's
   formula for lunar orientation (see [1]).
 
   In the case of DE-403, for the time period ~2000-2020, this
   time-dependent difference has an amplitude of approximately 0.005
   degrees, which is equivalent to approximately 150 m, measured along
   a great circle on the Moon's surface, while the average value is
   approximately 0.0025 degrees, or 75 m.


   Regarding Use of the ICRF in SPICE
   ==================================

   The IERS Celestial Reference Frame (ICRF) is offset from the J2000
   reference frame (equivalent to EME 2000) by a small rotation:  the
   J2000 pole offset magnitude is about 18 milliarcseconds (mas) and
   the equinox offset magnitude is approximately 78 milliarcseconds
   (see [3]).

   Certain SPICE data products *use the frame label "J2000" for data
   that actually are referenced to the ICRF.* This is the case for SPK
   files containing JPL version DE-4nn planetary ephemerides, for
   orientation data from generic text PCKs, and for binary PCKs,
   including binary lunar PCKs used in conjunction with this lunar frame
   kernel.

   Consequently, when SPICE computes the rotation between the "J2000"
   frame and the either of the lunar PA or ME frames, what's computed
   is actually the rotation between the ICRF and the respective
   lunar frame.

   Similarly, when a SPICE is used to compute the state given by a JPL
   DE planetary ephemeris SPK file of one ephemeris object relative to
   another (for example, the state of the Moon with respect to the
   Earth), expressed relative to the frame "J2000," the state is
   actually expressed relative to the ICRF.

   Because SPICE is already using the ICRF, users normally need not
   use the J2000-to-ICRF transformation to adjust results computed
   with SPICE.

   Lunar body-fixed frame associations
   =====================================================================

   By default, the SPICE system considers the body-fixed reference frame 
   associated with the Moon to be the one named IAU_MOON.  This body-frame
   association affects the outputs of the SPICE frame system routines

      CIDFRM
      CNMFRM
  
   and of the SPICE time conversion and geometry routines

      ET2LST
      ILLUM
      SRFXPT
      SUBPT
      SUBSOL

   Finally, any code that calls these routines to obtain results
   involving lunar body-fixed frames are affected.  Within SPICE, the
   only higher-level system that is affected is the dynamic frame
   system.

   NAIF provides "frame association" kernels that simplify changing the
   body-fixed frame associated with the Moon.  Using FURNSH to load
   either of the kernels named below sets the Moon's body-fixed frame 
   to that shown in the right-hand column:

      Kernel name          Lunar body-fixed frame        
      -----------          ----------------------
      moon_assoc_me.tf     MOON_ME
      moon_assoc_pa.tf     MOON_PA

   For further information see the in-line comments in the association
   kernels themselves.  Also see the Frames Required Reading section
   titled "Connecting an Object to its Body-fixed Frame."


   Using this Kernel
   =====================================================================

   In order for a SPICE-based application to use reference frames
   specified by this kernel, the application must load both this kernel
   and a binary lunar PCK containing lunar orientation data for the
   time of interest.  Normally the kernels need be loaded only once
   during program initialization.
 
   SPICE users may find it convenient to use a meta-kernel (also called
   a "FURNSH kernel") to name the kernels to be loaded. Below, we show
   an example of such a meta-kernel, as well as the source code of a
   small Fortran program that uses lunar body fixed frames.  The
   program's output is included as well.
 
   The kernel names shown here are simply used as examples; users must
   select the kernels appropriate for their applications.
 
   Numeric results shown below may differ very slightly from those obtained
   on users' computer systems.


   Meta-kernel
   -----------

      Example meta-kernel showing use of 

        - binary lunar PCK
        - generic lunar frame kernel (FK)
        - leapseconds kernel (LSK)
        - planetary SPK

       21-JUL-2006 (NJB)

       Note:  to actually use this kernel, replace the @ 
       characters below with backslashes (\).  The
       backslash character cannot be used here because these
       comments would be interpreted as actual load commands.

       This meta-kernel assumes that the referenced kernels exist
       in the user's current working directory.

          @begindata

            KERNELS_TO_LOAD = ( 'moon_pa_de403_1950-2198.bpc'
                                'moon_060721.tf'
                                'leapseconds.ker'
                                'de405.bsp'                  )

          @begintext


   Example program
   ---------------

            PROGRAM EX1
            IMPLICIT NONE

            INTEGER               FILSIZ
            PARAMETER           ( FILSIZ = 255 )

            CHARACTER*(FILSIZ)    META

            DOUBLE PRECISION      ET
            DOUBLE PRECISION      LT
            DOUBLE PRECISION      STME  ( 6 )
            DOUBLE PRECISION      STPA  ( 6 )

      C
      C     Prompt user for meta-kernel name.
      C
            CALL PROMPT ( 'Enter name of meta-kernel > ', META )

      C
      C     Load lunar PCK, generic lunar frame kernel,
      C     leapseconds kernel, and planetary ephemeris
      C     via metakernel.
      C
            CALL FURNSH ( META )

      C
      C     Convert a time of interest from UTC to ET.
      C
            CALL STR2ET ( '2006 jun 8 06:50:00', ET )

            WRITE (*,*) 'ET (sec past J2000 TDB): ', ET
            WRITE (*,*) '   State of Earth relative to Moon'

      C
      C     Find the geometric state of the Earth relative to the
      C     Moon at ET, expressed relative to the ME frame.
      C    
            CALL SPKEZR ( 'Earth',  ET,      'MOON_ME', 
           .              'NONE',   'Moon',  STME,      LT )

            WRITE (*,*) '      In MOON_ME frame:'
            WRITE (*,*) STME

      C
      C     Find the geometric state of the Earth relative to the
      C     Moon at ET, expressed relative to the PA frame.
      C    
            CALL SPKEZR ( 'Earth',  ET,      'MOON_PA', 
           .              'NONE',   'Moon',  STPA,      LT )

            WRITE (*,*) '      In MOON_PA frame:'
            WRITE (*,*) STPA

            END


   Program output
   --------------

   Enter name of meta-kernel > meta
    ET (sec past J2000 TDB):   203021465.
       State of Earth relative to Moon
          In MOON_ME frame:
     391739.183 -33210.254  25299.0887 -0.0592286405 -0.048721834  0.0917188552
          In MOON_PA frame:
     391719.148 -33331.588  25449.2934 -0.0592788895 -0.0487034073  0.0916961762



   References
   =====================================================================

   [1]  A.S. Konopliv, S.W. Asmar, E. Carranza, W.L. Sjogren, and D.N.
        Yuan (2001). "Recent Gravity Models as a Result of the Lunar
        Prospector Mission," Icarus 150, pp. 1-18.

   [2]  Seidelmann, P.K., Abalakin, V.K., Bursa, M., Davies, M.E.,
        Bergh, C. de, Lieske, J.H., Oberst, J., Simon, J.L., Standish,
        E.M., Stooke, P., and Thomas, P.C. (2002). "Report of the
        IAU/IAG Working Group on Cartographic Coordinates and Rotational
        Elements of the Planets and Satellites: 2000," Celestial
        Mechanics and Dynamical Astronomy, v.82, Issue 1, pp. 83-111.

   [3]  Roncoli, R. (2005).  "Lunar Constants and Models Document," 
        JPL D-32296.


   Frame Specifications
   =====================================================================

   MOON_PA is the name of the generic lunar principal axes (PA) reference
   frame.  This frame is an alias for the principal axes frame defined
   by the latest version of the JPL Solar System Dynamics Group's
   planetary ephemeris.
 
   In this instance of the lunar reference frames kernel, MOON_PA is an
   alias for the lunar principal axes frame associated with the
   planetary ephemeris DE-403.

   \begindata
 
      FRAME_MOON_PA                 = 31000
      FRAME_31000_NAME              = 'MOON_PA'
      FRAME_31000_CLASS             = 4
      FRAME_31000_CLASS_ID          = 31000
      FRAME_31000_CENTER            = 301

      TKFRAME_31000_SPEC            = 'MATRIX'
      TKFRAME_31000_RELATIVE        = 'MOON_PA_DE403'
      TKFRAME_31000_MATRIX          = ( 1 0 0
                                        0 1 0
                                        0 0 1 )

   \begintext

   MOON_ME is the name of the generic lunar mean Earth/ polar
   axis (ME) reference frame.  This frame is an alias for the mean
   Earth/polar axis frame defined by the latest version of the JPL
   Solar System Dynamics Group's planetary ephemeris.
 
   In this instance of the lunar reference frames kernel, MOON_ME is an
   alias for the lunar mean Earth/ polar axis frame associated with the
   planetary ephemeris DE-403.

   \begindata

      FRAME_MOON_ME                 = 31001
      FRAME_31001_NAME              = 'MOON_ME'
      FRAME_31001_CLASS             = 4
      FRAME_31001_CLASS_ID          = 31001
      FRAME_31001_CENTER            = 301
 
      TKFRAME_31001_SPEC            = 'MATRIX'
      TKFRAME_31001_RELATIVE        = 'MOON_ME_DE403'
      TKFRAME_31001_MATRIX          = ( 1 0 0
                                        0 1 0
                                        0 0 1 )


   \begintext


   MOON_PA_DE403 is the name of the lunar principal axes
   reference frame defined by JPL's DE-403 planetary ephemeris.

   \begindata

      FRAME_MOON_PA_DE403           = 31002
      FRAME_31002_NAME              = 'MOON_PA_DE403'
      FRAME_31002_CLASS             = 2
      FRAME_31002_CLASS_ID          = 31002
      FRAME_31002_CENTER            = 301


   \begintext

   MOON_ME_DE403 is the name of the lunar mean Earth/polar
   axis reference frame defined by JPL's DE-403 planetary ephemeris.

   Rotation angles are from reference [1].

   \begindata
 
      FRAME_MOON_ME_DE403           = 31003
      FRAME_31003_NAME              = 'MOON_ME_DE403'
      FRAME_31003_CLASS             = 4
      FRAME_31003_CLASS_ID          = 31003
      FRAME_31003_CENTER            = 301

      TKFRAME_31003_SPEC            = 'ANGLES'
      TKFRAME_31003_RELATIVE        = 'MOON_PA_DE403'
      TKFRAME_31003_ANGLES          = (   63.8986   79.0768   0.1462  )
      TKFRAME_31003_AXES            = (   3,        2,        1       )
      TKFRAME_31003_UNITS           = 'ARCSECONDS'

   \begintext


   Updating this Kernel
   --------------------

   When a new JPL DE providing lunar rotation data becomes available,
   the new lunar PA frame associated with that data set will be named

      MOON_PA_DEnnn

   where nnn is the version number of the DE.

   The PCK body ID code associated with that data set will be

      31004

   The frame ID and class ID for this frame will also be 31004.

   The generic PA frame specification will be updated to point to the
   new DE-specific PA frame.  The rest of this frame specification
   is unchanged.

   The ME frame name associated with the new data set will be named
 
      MOON_ME_DEnnn

   The frame ID and class ID for this frame will be 

      31005

   The rotational offset between this frame and the new DE-specific PA
   frame will need to be updated; this offset is DE-dependent.
    
   The generic ME frame specification will be updated to point to the
   new DE-specific ME frame.  The rest of this frame specification
   is unchanged.



   =====================================================================
   End of kernel









