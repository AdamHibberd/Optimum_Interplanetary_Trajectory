
KPL/FK


   Topocentric Reference Frame Definition Kernel for DSN Stations
   =====================================================================

   Original file name:                   earth_topo_050714.tf
   Creation date:                        2005 July 14 21:00
   Created by:                           Nat Bachman  (NAIF/JPL)


   Introduction
   =====================================================================

   This file defines topocentric reference frames associated with each
   of the DSN stations cited in the list below under "Position Data."
   Each topocentric reference frame ("frame" for short) is centered at
   the associated station and is fixed to the earth. Mathematically, a
   frame "definition" is a specification of the orientation of the
   frame relative to another frame.  In this file, the other frame,
   which we'll refer to as the "base frame," is the terrestrial
   reference frame ITRF93.
 
   The orientation of a topocentric frame relative to the base frame
   relies on a reference spheroid (see "Data Sources" below).  The
   z-axis of the topocentric frame contains the station location and is
   normal to the reference spheroid:  the line containing the z-axis
   intersects the reference spheroid at right angles. The x-axis
   points north and the y-axis points west.  Note that stations
   normally have non-zero altitude with respect to the spheroid.

   Loosely speaking, a topocentric frame enables computations involving
   the local directions "north", "west," and "up" at a surface point on
   an extended body. For example, the "elevation" of an object relative
   to the center of a topocentric frame is the object's colatitude in
   that frame.  The corresponding azimuth is the angle from the
   topocentric frame's x-axis to the projection of the center-to-object
   vector onto the topocentric frame's x-y plane, measured in the
   clockwise direction.

   The orientation of a topocentric frame relative to the base frame can
   be described by an Euler angle sequence.  Let M be the rotation
   matrix that maps vectors from the base frame to a specified
   topocentric frame.  Then

                
      M   =  [ Pi  ]  [ Pi/2 - LAT ]  [ LON ]
                    3               2        3

   where LON, LAT are the associated station's geodetic latitude and
   longitude.  Note that the frame definitions below actually 
   provide Euler angles for the inverse of M and use units of
   degrees, so the angle sequences are

       -1                         o          o
      M   =  [ -LON ]   [ LAT - 90 ]    [ 180 ]
                     3               2         3

   See the Rotation Required Reading for details concerning Euler
   angles.


   Using this kernel
   =====================================================================

   Revision description
   --------------------

   This kernel supersedes  

      earth_topo_040916.tf

   This revision defines an additional topocentric reference frame
   centered at the new position of the relocated DSS-65 antenna; the
   antenna has been moved by about 61m.

   Some users of station location data are adopting the name

      DSS-64 

   to represent the new location; others are continuing to use
   the old name

      DSS-65

   This kernel enables SPICE users to refer to a topocentric frame
   centered at the new location by either of the names

      DSS-65_TOPO
      DSS-64_TOPO

   This kernel does not define a topocentric frame centered at the
   old location of DSS-65:  the previous version of this kernel
   may be used to provide that definition.
 
 
   Planned updates
   ---------------

   NAIF plans to replace this kernel with one containing additional
   data for tracking stations at Noto and New Norcia.  Data for the
   sites covered by this file will be unchanged in that update.


   Kernel loading
   --------------

   In order for a SPICE-based program to make use of this kernel,
   the kernel must be loaded via the SPICE routine FURNSH.  If you
   are running application software created by a third party, see the
   documentation for that software for instructions on kernel
   management.

   See also "Associated SPK files" and "Associated PCK files" below.


   Base frame alias
   ----------------

   This kernel uses the frame alias 'EARTH_FIXED' to designate the
   base frame.  Below, this alias is mapped to the frame name 'ITRF93'.
   In some situations, for example when low accuracy, long term 
   predictions are desired, it may be convenient to map EARTH_FIXED
   to 'IAU_EARTH'.

   See the Frames Required Reading for details.


   Associated PCK files
   --------------------

   For high-accuracy work, this kernel should be used together with a
   high-precision, binary earth PCK file.  

      NAIF produces these kernels on a regular basis; they can be
      obtained via anonymous ftp from the NAIF server

         naif.jpl.nasa.gov

      The PCK is located in the path

         pub/naif/generic_kernels/pck

      The file name is of the form

         earth_000101_yymmdd_yymmdd.bpc

      The first two dates are the file's start and stop times; the third
      is the epoch of the last datum in the EOP file:  data from
      this epoch forward are predicted.

      The file's coverage starts at a fixed date (currently chosen to
      be 2000 Jan. 1) and extends to the end of the predict region,
      which has a duration of roughly 3 months.


   For less accurate work, a text PCK may suffice.  To use this kernel
   with a text PCK, the base frame alias EARTH_FIXED must be mapped to
   'IAU_EARTH'.  Text PCKs may be appropriate for work involving
   long term predicts.


   Associated SPK files
   --------------------

   This file is compatible with the SPK files
 
       earthstns_fx_050714.bsp       [reference frame: EARTH_FIXED]
       earthstns_itrf93_050714.bsp   [reference frame: ITRF93     ]

   both of which provide state vectors for each station covered by this
   file.

   Most applications will need to load one of the above SPK files in
   order to make use of this frame kernel.


   DSS-64 and DSS-65
   -----------------

   See "Revision description" above for a description of the data 
   coverage provided by this file for DSS-64 and DSS-65.

   To enable use of the name DSS-64, user applications may load
   a text kernel containing the assignments
 
 
          NAIF_BODY_NAME  +=  'DSS-64'
          NAIF_BODY_CODE  +=   399064
 
 
   This frame kernel includes the necessary definitions.
 
   See the NAIF_IDs Required Reading for details.
 
 
   PARKES
   ------
 
   The station location data source produced by JPL's section 335
   now refers to the Parkes station as "DSS-49."  The SPICE Toolkit
   currently supports the NAIF ID code/name mappings
 
      399005  <-->  DSS-05  (secondary)
              <-->  PARKES  (primary)
 
   Identical ephemeris data are provided in this file for both ID codes
   399005 and 399049.
 
   See the NAIF_IDs Required Reading for details.
 
   In this file, all of the frame names

      DSS-05_TOPO
      PARKES_TOPO
      DSS-49_TOPO

   are associated with the topocentric frame for Parkes, so any of
   these names will be recognized by the SPICE Toolkit when this kernel
   is loaded.  All three names refer to mathematically equivalent
   frames. 



   Data sources
   =====================================================================

   The data described here are taken from the JPL web site at URL

      http://epic/nav/eop/stations.html

   The site is maintained by Tod Ratcliff, JPL section 335.

   Additional source:

      Location data for DSS-64 are from an e-mail communication 
      from W. M. Folkner to N. J. Bachman, dated June 23, 2005. 


   Reference Spheroid
   ------------------

   The reference bi-axial spheroid is defined by an equatorial and a
   polar radius.  Calling these Re and Rp respectively, the flattening
   factor f is defined as
  
      f = ( Re - Rp ) / Re
   
   For the reference spheroid used by this file, the equatorial radius
   Re and inverse flattening factor 1/f are

      Re  = 6378136.3 m
      1/f = 298.257


   Position data
   -------------

   The Cartesian station locations from which the topocentric frames
   defined here were derived are shown below. Station locations in the
   ITRF93 frame are:

       Antenna  Diameter   x (m)            y (m)           z (m)
       DSS 12   34m   -2350444.0057   -4651980.7620    3665630.9322
       DSS 13   34m   -2351112.6586   -4655530.6359    3660912.7276
       DSS 14   70m   -2353621.4197   -4641341.4717    3677052.3178
       DSS 15   34m   -2353538.9575   -4641649.4287    3676669.9837
       DSS 16   26m   -2354763.3257   -4646787.3837    3669387.0099
       DSS 17    9m   -2354730.5247   -4646751.6975    3669440.5998
       DSS 23   11m   -2354757.7341   -4646934.5965    3669207.7651
       DSS 24   34m   -2354906.7087   -4646840.0834    3669242.3207
       DSS 25   34m   -2355022.0140   -4646953.2040    3669040.5666
       DSS 26   34m   -2354890.7996   -4647166.3182    3668871.7546
       DSS 27   34m   -2349915.4275   -4656756.4059    3660096.4693
       DSS 28   34m   -2350102.0169   -4656673.3686    3660103.5180
       DSS 33   11m   -4461083.8425    2682281.6961   -3674569.9725
       DSS 34   34m   -4461147.0925    2682439.2385   -3674393.1332
       DSS 42   34m   -4460981.3463    2682413.4680   -3674581.6534
       DSS 43   70m   -4460894.9170    2682361.5070   -3674748.1517
       DSS 45   34m   -4460935.5783    2682765.6611   -3674380.9824
       DSS 46   26m   -4460828.9473    2682129.5071   -3674975.0884
       DSS 49   64m   -4554232.1933    2816758.9161   -3454035.6434
       DSS 53   11m    4849330.0161    -360337.8678    4114758.9123
       DSS 54   34m    4849434.4877    -360723.8999    4114618.8354
       DSS 55   34m    4849525.2561    -360606.0932    4114495.0843
       DSS 61   34m    4849245.0787    -360277.9478    4114884.5772
       DSS 63   70m    4849092.5175    -360180.3480    4115109.2506
       DSS 64   34m    4849339.6448    -360427.6560    4114750.7428

       DSS 65   34m    Prior to July 3, 2005:
                       4849336.6176    -360488.6349    4114748.9218

                       After July 3, 2005:
                       4849339.6448    -360427.6560    4114750.7428
                       (Same as DSS 64)

       DSS 66   26m    4849148.4311    -360474.6175    4114995.1679


   Epoch
   -----

   The epoch associated with these data is given by the source as
   "2003.0."  The time variation of the data is slow enough so that
   specification of the time system is unimportant. However, in the
   creation of this file, the epoch is assumed to be

      2003 Jan 1 00:00:00 TDB

   The movement of the stations due to tectonic plate motion is taken
   into account in creation of the frame definitions used in this file:
   the center locations and orientations of the reference frames are 
   associated with station locations extrapolated to the date

      2005 July 15 00:00:00 TDB

   This extrapolation results in a small rotation of the frames relative
   to their orientations as given by the previous version of this kernel.
   The rotations are typically on the order of 20 nanoradians.


   Accuracy
   --------

   The Euler angles specified in this kernel have much lower
   accuracy than suggested by their presentation as double-precision
   numbers.  The presentation was selected to avoid loss of precision.

   The frame definitions given here correspond to station locations
   at a fixed epoch.  Because station locations are time-varying,
   this kernel will gradually become inconsistent with the 
   corresponding station location data.

   The following discussion concerning station location accuracy is
   from the referenced web site.  The citations in the text refer to
   the documents listed below.

      The uncertainty in the station locations is described by a
      covariance matrix [2]. The coordinate uncertainties, given by the
      square-root of the diagonal elements of the covariance matrix,
      are about 4 cm for DSN stations which have participated in
      regular VLBI experiments, and about 10 cm for other stations.
      These coordinate uncertaintes [sic] do not account for
      uncertainties in Earth orientation calibrations. Uncertainties in
      Earth orientation as applied to spacecraft navigation are
      discussed in [5].


   References
   ----------

   The site lists the following references:

       1. C. Boucher, Z. Altamimi, L. Duhem, 
          "Results and analysis of the ITRF93", 
          IERS Technical Note 18, Observatoire de Paris, 1994.

       2. W. M. Folkner, DSN station locations and uncertainties, 
          JPL TDA Progress Report, 42-128,pp. 1-34,1996.

       3. T. Moyer, "Mathematical formulation of the double-precision 
          Orbit Determination Program", JPL Technical Report 32-1527, 1971

       4. C. S. Jacobs and A. Rius, Internal consistency of VLBI surveying 
          between DSS 63 and DSS 65", JPL IOM 335.3-90-034, 11 May 1992.

       5. J. A. Estefan and W. M. Folkner, Sensitivity of planetary cruise 
          navigation to Earth orientation calibration errors, JPL TDA 
          Progress Report 42-123, pp. 1-29, 1995.

       6. T. D. Moyer, "Frame tie rotations and nutation corrections for 
          the ODP", JPL EM 314-558, 26 February 1993.

       7. E. M. Standish, X X Newhall, J. G. Williams, W. M. Folkner, 
          "JPL planetary and lunar ephemerides DE403/LE403", JPL IOM 
          314.10-127, 22 May 1995.




   Reference frame definitions
   =====================================================================


   EARTH_FIXED alias mapping
   -------------------------

   Constant-offset frame definition for the frame alias EARTH_FIXED:
   EARTH_FIXED is mapped to ITRF93.


\begindata

   TKFRAME_EARTH_FIXED_RELATIVE = 'ITRF93'
   TKFRAME_EARTH_FIXED_SPEC     = 'MATRIX'
   TKFRAME_EARTH_FIXED_MATRIX   = ( 1   0   0
                                    0   1   0
                                    0   0   1 )

\begintext

   Name-ID code associations
   --------------------------------

   PARKES is called "DSS-49" in the data source.  The older names DSS-05
   and PARKES are associated with the ID code 399005 for backward
   compatibility.
 
   The ID code 399064 is not yet a SPICE built-in code, so it is
   associated here with the name DSS-64.

\begindata

   NAIF_BODY_NAME                        +=  'DSS-64'
   NAIF_BODY_CODE                        +=  399064

\begintext

   DSN frame definitions follow.

\begindata

   FRAME_PARKES_TOPO                     = 1399005
   FRAME_1399005_NAME                    = 'PARKES_TOPO'
   FRAME_1399005_CLASS                   = 4
   FRAME_1399005_CLASS_ID                = 1399005
   FRAME_1399005_CENTER                  = 399005

   OBJECT_399005_FRAME                   = 'PARKES_TOPO'

   TKFRAME_PARKES_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_PARKES_TOPO_SPEC              = 'ANGLES'
   TKFRAME_PARKES_TOPO_UNITS             = 'DEGREES'
   TKFRAME_PARKES_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_PARKES_TOPO_ANGLES            = ( -148.2635161528947,
                                             -122.9983980423326,
                                              180.0000000000000  )


   FRAME_DSS-12_TOPO                     = 1399012
   FRAME_1399012_NAME                    = 'DSS-12_TOPO'
   FRAME_1399012_CLASS                   = 4
   FRAME_1399012_CLASS_ID                = 1399012
   FRAME_1399012_CENTER                  = 399012

   OBJECT_399012_FRAME                   = 'DSS-12_TOPO'

   TKFRAME_DSS-12_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-12_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-12_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-12_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-12_TOPO_ANGLES            = ( -243.1945102442646,
                                              -54.7000629043147,
                                              180.0000000000000  )


   FRAME_DSS-13_TOPO                     = 1399013
   FRAME_1399013_NAME                    = 'DSS-13_TOPO'
   FRAME_1399013_CLASS                   = 4
   FRAME_1399013_CLASS_ID                = 1399013
   FRAME_1399013_CENTER                  = 399013

   OBJECT_399013_FRAME                   = 'DSS-13_TOPO'

   TKFRAME_DSS-13_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-13_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-13_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-13_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-13_TOPO_ANGLES            = ( -243.2055404763616,
                                              -54.7528357325366,
                                              180.0000000000000  )


   FRAME_DSS-14_TOPO                     = 1399014
   FRAME_1399014_NAME                    = 'DSS-14_TOPO'
   FRAME_1399014_CLASS                   = 4
   FRAME_1399014_CLASS_ID                = 1399014
   FRAME_1399014_CENTER                  = 399014

   OBJECT_399014_FRAME                   = 'DSS-14_TOPO'

   TKFRAME_DSS-14_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-14_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-14_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-14_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-14_TOPO_ANGLES            = ( -243.1104612607222,
                                              -54.5740991182250,
                                              180.0000000000000  )


   FRAME_DSS-15_TOPO                     = 1399015
   FRAME_1399015_NAME                    = 'DSS-15_TOPO'
   FRAME_1399015_CLASS                   = 4
   FRAME_1399015_CLASS_ID                = 1399015
   FRAME_1399015_CENTER                  = 399015

   OBJECT_399015_FRAME                   = 'DSS-15_TOPO'

   TKFRAME_DSS-15_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-15_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-15_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-15_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-15_TOPO_ANGLES            = ( -243.1128043663782,
                                              -54.5781467088189,
                                              180.0000000000000  )


   FRAME_DSS-16_TOPO                     = 1399016
   FRAME_1399016_NAME                    = 'DSS-16_TOPO'
   FRAME_1399016_CLASS                   = 4
   FRAME_1399016_CLASS_ID                = 1399016
   FRAME_1399016_CENTER                  = 399016

   OBJECT_399016_FRAME                   = 'DSS-16_TOPO'

   TKFRAME_DSS-16_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-16_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-16_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-16_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-16_TOPO_ANGLES            = ( -243.1263497222477,
                                              -54.6584605941241,
                                              180.0000000000000  )


   FRAME_DSS-17_TOPO                     = 1399017
   FRAME_1399017_NAME                    = 'DSS-17_TOPO'
   FRAME_1399017_CLASS                   = 4
   FRAME_1399017_CLASS_ID                = 1399017
   FRAME_1399017_CENTER                  = 399017

   OBJECT_399017_FRAME                   = 'DSS-17_TOPO'

   TKFRAME_DSS-17_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-17_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-17_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-17_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-17_TOPO_ANGLES            = ( -243.1264941091303,
                                              -54.6578234080814,
                                              180.0000000000000  )


   FRAME_DSS-23_TOPO                     = 1399023
   FRAME_1399023_NAME                    = 'DSS-23_TOPO'
   FRAME_1399023_CLASS                   = 4
   FRAME_1399023_CLASS_ID                = 1399023
   FRAME_1399023_CENTER                  = 399023

   OBJECT_399023_FRAME                   = 'DSS-23_TOPO'

   TKFRAME_DSS-23_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-23_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-23_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-23_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-23_TOPO_ANGLES            = ( -243.1271364494682,
                                              -54.6604496329255,
                                              180.0000000000000  )


   FRAME_DSS-24_TOPO                     = 1399024
   FRAME_1399024_NAME                    = 'DSS-24_TOPO'
   FRAME_1399024_CLASS                   = 4
   FRAME_1399024_CLASS_ID                = 1399024
   FRAME_1399024_CENTER                  = 399024

   OBJECT_399024_FRAME                   = 'DSS-24_TOPO'

   TKFRAME_DSS-24_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-24_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-24_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-24_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-24_TOPO_ANGLES            = ( -243.1252050470835,
                                              -54.6601071639813,
                                              180.0000000000000  )


   FRAME_DSS-25_TOPO                     = 1399025
   FRAME_1399025_NAME                    = 'DSS-25_TOPO'
   FRAME_1399025_CLASS                   = 4
   FRAME_1399025_CLASS_ID                = 1399025
   FRAME_1399025_CENTER                  = 399025

   OBJECT_399025_FRAME                   = 'DSS-25_TOPO'

   TKFRAME_DSS-25_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-25_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-25_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-25_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-25_TOPO_ANGLES            = ( -243.1246362656502,
                                              -54.6623880235187,
                                              180.0000000000000  )


   FRAME_DSS-26_TOPO                     = 1399026
   FRAME_1399026_NAME                    = 'DSS-26_TOPO'
   FRAME_1399026_CLASS                   = 4
   FRAME_1399026_CLASS_ID                = 1399026
   FRAME_1399026_CENTER                  = 399026

   OBJECT_399026_FRAME                   = 'DSS-26_TOPO'

   TKFRAME_DSS-26_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-26_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-26_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-26_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-26_TOPO_ANGLES            = ( -243.1269829764524,
                                              -54.6643107688601,
                                              180.0000000000000  )


   FRAME_DSS-27_TOPO                     = 1399027
   FRAME_1399027_NAME                    = 'DSS-27_TOPO'
   FRAME_1399027_CLASS                   = 4
   FRAME_1399027_CLASS_ID                = 1399027
   FRAME_1399027_CENTER                  = 399027

   OBJECT_399027_FRAME                   = 'DSS-27_TOPO'

   TKFRAME_DSS-27_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-27_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-27_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-27_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-27_TOPO_ANGLES            = ( -243.2233490223594,
                                              -54.7617282056820,
                                              180.0000000000000  )


   FRAME_DSS-28_TOPO                     = 1399028
   FRAME_1399028_NAME                    = 'DSS-28_TOPO'
   FRAME_1399028_CLASS                   = 4
   FRAME_1399028_CLASS_ID                = 1399028
   FRAME_1399028_CENTER                  = 399028

   OBJECT_399028_FRAME                   = 'DSS-28_TOPO'

   TKFRAME_DSS-28_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-28_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-28_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-28_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-28_TOPO_ANGLES            = ( -243.2211082995305,
                                              -54.7617279659690,
                                              180.0000000000000  )


   FRAME_DSS-33_TOPO                     = 1399033
   FRAME_1399033_NAME                    = 'DSS-33_TOPO'
   FRAME_1399033_CLASS                   = 4
   FRAME_1399033_CLASS_ID                = 1399033
   FRAME_1399033_CENTER                  = 399033

   OBJECT_399033_FRAME                   = 'DSS-33_TOPO'

   TKFRAME_DSS-33_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-33_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-33_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-33_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-33_TOPO_ANGLES            = ( -148.9830923595113,
                                             -125.4004837705609,
                                              180.0000000000000  )


   FRAME_DSS-34_TOPO                     = 1399034
   FRAME_1399034_NAME                    = 'DSS-34_TOPO'
   FRAME_1399034_CLASS                   = 4
   FRAME_1399034_CLASS_ID                = 1399034
   FRAME_1399034_CENTER                  = 399034

   OBJECT_399034_FRAME                   = 'DSS-34_TOPO'

   TKFRAME_DSS-34_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-34_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-34_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-34_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-34_TOPO_ANGLES            = ( -148.9819650021110,
                                             -125.3984778756552,
                                              180.0000000000000  )


   FRAME_DSS-42_TOPO                     = 1399042
   FRAME_1399042_NAME                    = 'DSS-42_TOPO'
   FRAME_1399042_CLASS                   = 4
   FRAME_1399042_CLASS_ID                = 1399042
   FRAME_1399042_CENTER                  = 399042

   OBJECT_399042_FRAME                   = 'DSS-42_TOPO'

   TKFRAME_DSS-42_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-42_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-42_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-42_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-42_TOPO_ANGLES            = ( -148.9812679796550,
                                             -125.4006736665826,
                                              180.0000000000000  )


   FRAME_DSS-43_TOPO                     = 1399043
   FRAME_1399043_NAME                    = 'DSS-43_TOPO'
   FRAME_1399043_CLASS                   = 4
   FRAME_1399043_CLASS_ID                = 1399043
   FRAME_1399043_CENTER                  = 399043

   OBJECT_399043_FRAME                   = 'DSS-43_TOPO'

   TKFRAME_DSS-43_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-43_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-43_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-43_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-43_TOPO_ANGLES            = ( -148.9812678907116,
                                             -125.4024232666378,
                                              180.0000000000000  )


   FRAME_DSS-45_TOPO                     = 1399045
   FRAME_1399045_NAME                    = 'DSS-45_TOPO'
   FRAME_1399045_CLASS                   = 4
   FRAME_1399045_CLASS_ID                = 1399045
   FRAME_1399045_CENTER                  = 399045

   OBJECT_399045_FRAME                   = 'DSS-45_TOPO'

   TKFRAME_DSS-45_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-45_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-45_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-45_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-45_TOPO_ANGLES            = ( -148.9776862148021,
                                             -125.3984567187688,
                                              180.0000000000000  )


   FRAME_DSS-46_TOPO                     = 1399046
   FRAME_1399046_NAME                    = 'DSS-46_TOPO'
   FRAME_1399046_CLASS                   = 4
   FRAME_1399046_CLASS_ID                = 1399046
   FRAME_1399046_CENTER                  = 399046

   OBJECT_399046_FRAME                   = 'DSS-46_TOPO'

   TKFRAME_DSS-46_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-46_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-46_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-46_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-46_TOPO_ANGLES            = ( -148.9830822665351,
                                             -125.4050096708302,
                                              180.0000000000000  )


   FRAME_DSS-49_TOPO                     = 1399049
   FRAME_1399049_NAME                    = 'DSS-49_TOPO'
   FRAME_1399049_CLASS                   = 4
   FRAME_1399049_CLASS_ID                = 1399049
   FRAME_1399049_CENTER                  = 399049

   OBJECT_399049_FRAME                   = 'DSS-49_TOPO'

   TKFRAME_DSS-49_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-49_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-49_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-49_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-49_TOPO_ANGLES            = ( -148.2635161528947,
                                             -122.9983980423326,
                                              180.0000000000000  )


   FRAME_DSS-53_TOPO                     = 1399053
   FRAME_1399053_NAME                    = 'DSS-53_TOPO'
   FRAME_1399053_CLASS                   = 4
   FRAME_1399053_CLASS_ID                = 1399053
   FRAME_1399053_CENTER                  = 399053

   OBJECT_399053_FRAME                   = 'DSS-53_TOPO'

   TKFRAME_DSS-53_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-53_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-53_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-53_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-53_TOPO_ANGLES            = ( -355.7503485315841,
                                              -49.5726421412624,
                                              180.0000000000000  )


   FRAME_DSS-54_TOPO                     = 1399054
   FRAME_1399054_NAME                    = 'DSS-54_TOPO'
   FRAME_1399054_CLASS                   = 4
   FRAME_1399054_CLASS_ID                = 1399054
   FRAME_1399054_CENTER                  = 399054

   OBJECT_399054_FRAME                   = 'DSS-54_TOPO'

   TKFRAME_DSS-54_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-54_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-54_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-54_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-54_TOPO_ANGLES            = ( -355.7459038709393,
                                              -49.5743777507050,
                                              180.0000000000000  )


   FRAME_DSS-55_TOPO                     = 1399055
   FRAME_1399055_NAME                    = 'DSS-55_TOPO'
   FRAME_1399055_CLASS                   = 4
   FRAME_1399055_CLASS_ID                = 1399055
   FRAME_1399055_CENTER                  = 399055

   OBJECT_399055_FRAME                   = 'DSS-55_TOPO'

   TKFRAME_DSS-55_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-55_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-55_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-55_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-55_TOPO_ANGLES            = ( -355.7473673994048,
                                              -49.5757035288805,
                                              180.0000000000000  )


   FRAME_DSS-61_TOPO                     = 1399061
   FRAME_1399061_NAME                    = 'DSS-61_TOPO'
   FRAME_1399061_CLASS                   = 4
   FRAME_1399061_CLASS_ID                = 1399061
   FRAME_1399061_CENTER                  = 399061

   OBJECT_399061_FRAME                   = 'DSS-61_TOPO'

   TKFRAME_DSS-61_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-61_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-61_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-61_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-61_TOPO_ANGLES            = ( -355.7509784608760,
                                              -49.5712602628462,
                                              180.0000000000000  )


   FRAME_DSS-63_TOPO                     = 1399063
   FRAME_1399063_NAME                    = 'DSS-63_TOPO'
   FRAME_1399063_CLASS                   = 4
   FRAME_1399063_CLASS_ID                = 1399063
   FRAME_1399063_CENTER                  = 399063

   OBJECT_399063_FRAME                   = 'DSS-63_TOPO'

   TKFRAME_DSS-63_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-63_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-63_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-63_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-63_TOPO_ANGLES            = ( -355.7519921564240,
                                              -49.5687896830942,
                                              180.0000000000000  )


   FRAME_DSS-64_TOPO                     = 1399064
   FRAME_1399064_NAME                    = 'DSS-64_TOPO'
   FRAME_1399064_CLASS                   = 4
   FRAME_1399064_CLASS_ID                = 1399064
   FRAME_1399064_CENTER                  = 399064

   OBJECT_399064_FRAME                   = 'DSS-64_TOPO'

   TKFRAME_DSS-64_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-64_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-64_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-64_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-64_TOPO_ANGLES            = ( -355.7493019024685,
                                              -49.5727930629388,
                                              180.0000000000000  )


   FRAME_DSS-65_TOPO                     = 1399065
   FRAME_1399065_NAME                    = 'DSS-65_TOPO'
   FRAME_1399065_CLASS                   = 4
   FRAME_1399065_CLASS_ID                = 1399065
   FRAME_1399065_CENTER                  = 399065

   OBJECT_399065_FRAME                   = 'DSS-65_TOPO'

   TKFRAME_DSS-65_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-65_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-65_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-65_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-65_TOPO_ANGLES            = ( -355.7493019024685,
                                              -49.5727930629388,
                                              180.0000000000000  )


   FRAME_DSS-66_TOPO                     = 1399066
   FRAME_1399066_NAME                    = 'DSS-66_TOPO'
   FRAME_1399066_CLASS                   = 4
   FRAME_1399066_CLASS_ID                = 1399066
   FRAME_1399066_CENTER                  = 399066

   OBJECT_399066_FRAME                   = 'DSS-66_TOPO'

   TKFRAME_DSS-66_TOPO_RELATIVE          = 'EARTH_FIXED'
   TKFRAME_DSS-66_TOPO_SPEC              = 'ANGLES'
   TKFRAME_DSS-66_TOPO_UNITS             = 'DEGREES'
   TKFRAME_DSS-66_TOPO_AXES              = ( 3, 2, 3)
   TKFRAME_DSS-66_TOPO_ANGLES            = ( -355.7485830705283,
                                              -49.5700245583355,
                                              180.0000000000000  )

\begintext

