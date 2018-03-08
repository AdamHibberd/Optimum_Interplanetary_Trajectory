KPL/FK

\beginlabel
PDS_VERSION_ID               = PDS3
RECORD_TYPE                  = STREAM
RECORD_BYTES                 = "N/A"
^SPICE_KERNEL                = "di_v17.tf"
MISSION_NAME                 = "DEEP IMPACT"
SPACECRAFT_NAME              = {
                               "DEEP IMPACT FLYBY SPACECRAFT",
                               "DEEP IMPACT IMPACTOR SPACECRAFT"
                               }
DATA_SET_ID                  = "DI-C-SPICE-6-V1.0"
KERNEL_TYPE_ID               = FK
PRODUCT_ID                   = "di_v17.tf"
PRODUCT_CREATION_TIME        = 2006-12-21T10:53:28
PRODUCER_ID                  = "NAIF/JPL"
MISSION_PHASE_NAME           = "N/A"
PRODUCT_VERSION_TYPE         = ACTUAL
PLATFORM_OR_MOUNTING_NAME    = "N/A"
START_TIME                   = "N/A"
STOP_TIME                    = "N/A"
SPACECRAFT_CLOCK_START_COUNT = "N/A"
SPACECRAFT_CLOCK_STOP_COUNT  = "N/A"
TARGET_NAME                  = "9P/TEMPEL 1 (1867 G1)"
INSTRUMENT_NAME              = "N/A"
NAIF_INSTRUMENT_ID           = "N/A"
SOURCE_PRODUCT_ID            = "N/A"
NOTE                         = "See comments in the file for details"
OBJECT                       = SPICE_KERNEL
  INTERCHANGE_FORMAT         = ASCII
  KERNEL_TYPE                = FRAMES
  DESCRIPTION                = "DI SPICE Frames Kernel file providing frame
definitions for spacecraft and instrument frames for both spacecraft (DIF
and DII) as well as the body-fixed frame definition for 9P/Tempel 1, created
by NAIF, JPL. "
END_OBJECT                   = SPICE_KERNEL
\endlabel


Deep Impact Flyby and Impactor Spacecraft Frames Kernel
===============================================================================

   This frame kernel contains complete set of frame definitions for the 
   Deep Impact Flyby (DIF) and Deep Impact Impactor (DII) spacecraft 
   including definitions for the s/c fixed frames and science instrument 
   frames. This kernel also contains NAIF ID/name mapping for the Deep
   Impact spacecraft and instruments.


Version and Date
-------------------------------------------------------------------------------

   Version 1.7 -- August 1, 2006 -- Boris Semenov, NAIF

      Incorporated HRI/IR misalignment from [12].

   Version 1.6 -- September 12, 2005 -- Boris Semenov, NAIF

      Added DIF_HGA_GIMBAL frame (frame ID -140930). Changed DIF_HGA
      frame ID to -140930 and redefined it as a fixed offset frame 
      w.r.t. to the DIF_HGA_GIMBAL.

   Version 1.5 -- August 1, 2005 -- Boris Semenov, NAIF

      Corrected typo in the DII_ITS_ADCS definition. The CENTER was 
      incorrectly set to -140. It was changed to -70.

   Version 1.4 -- June 28, 2005 -- Boris Semenov, NAIF

      Incorporated final MRI misalignment values (see [11]).

   Version 1.3 -- April 5, 2005 -- Boris Semenov, NAIF

      Fixed typo in the comments in ITS section.

   Version 1.2 -- February 26, 2005 -- Boris Semenov, NAIF

      For completeness, defined additional HRI, MRI, and ITS frames
      following ADCS convention. All three -- DIF_HRI_VIS_ADCS,
      DIF_MRI_ADCS, and DII_ITS_ADCS -- are defined as fixed offset,
      nominally-rotated frames relative to the corresponding camera
      frames.

   Version 1.1 -- February 11, 2005 -- Boris Semenov, NAIF

      Incorporated initial OPNAV estimate of the MRI misalignment.

   Version 1.0 -- December 20, 2004 -- Boris Semenov, NAIF

      Re-defined HRI, MRI, and ITS frames to be consistent with OPNAV
      and science expectations. Described impactor and ITS frames.

   Version 0.3 -- April 28, 2003 -- Boris Semenov, NAIF

      Added HGA frame.

   Version 0.2 -- September 17, 2002 -- Boris Semenov, NAIF

      Added LGA frames (preliminary; antenna +X is not aligned with the
      pattern reference direction.)

   Version 0.1 -- April 16, 2002 -- Boris Semenov, NAIF

      Added star tracker (ST#) and coarse solar sensor (CSS#) frames.

   Version 0.0 -- April 13, 2001 -- Boris Semenov, NAIF

      Preliminary Version: 

           ---   TO BE USED ONLY FOR SOA TOOL DEMONSTRATION    --- 

      This version is based solely on [4], which, as of Rev 010216, did 
      NOT contain complete set of information required to understand and
      define Deep Impact Flyby and Impactor spacecraft and their 
      instrument frames without making certain "guesses."  
      

References
-------------------------------------------------------------------------------

   1. ``Frames Required Reading''

   2. ``Kernel Pool Required Reading''

   3. ``C-Kernel Required Reading''

   4. ``Deep Impact Instruments Requirements Specification'', 2001
      February 16, M. Ensminger

   5. E-mail from L. Kendall to M. Hughes providing star tracker and 
      coarse sun sensor boresight directions, March 27, 2002

   6. ``Deep Impact Coarse Sun Sensor Placement'', Systems Engineering 
      Report, by M. Larson, April 2002 (includes changes per ECR from 
      December 17, 2001)

   7. ``Deep Impact LGA Pattern Analysis'', DI-SC-COM-063, 07/17/02

   8. ``Deep Impact Camera Orientations: AutoNav and ADCS
      Definitions'', March 10, 2003

   9. E-mail from Brian Carcich, December 20, 2004.
   
   10. E-mail from Nick Mastrodemos, DI OPNAV, JPL, January 18, 2005.

   11. "Deep Impact Precise Alignments Summary", Rev. A: 09 June 2005, by
      Steve Waydo, DI AACS.

   12. DI Calibration Document (DICalPaper062706_inst_alignments.pdf)


Contact Information
-------------------------------------------------------------------------------

   Boris V. Semenov, NAIF/JPL, (818)-354-8136, bsemenov@spice.jpl.nasa.gov


Implementation Notes
-------------------------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make
   use of this frame kernel must `load' the kernel, normally during program
   initialization. The SPICELIB routine FURNSH loads a kernel file into the 
   pool as shown below.

      CALL FURNSH ( frame_kernel_name )

   This file was created and may be updated with a text editor or word
   processor.


Deep Impact NAIF ID Codes
-------------------------------------------------------------------------------

   Deep Impact Flyby (DIF) spacecraft and instruments IDs:

   \begindata

      NAIF_BODY_NAME += ( 'DEEP_IMPACT_FLYBY_SC'    )
      NAIF_BODY_CODE += ( -140                      )

      NAIF_BODY_NAME += ( 'DIF'                     )
      NAIF_BODY_CODE += ( -140                      )

      NAIF_BODY_NAME += ( 'DIF_HRI'                 )
      NAIF_BODY_CODE += ( -140100                   )

      NAIF_BODY_NAME += ( 'DIF_HRI_VIS'             )
      NAIF_BODY_CODE += ( -140110                   )

      NAIF_BODY_NAME += ( 'DIF_HRI_IR'              )
      NAIF_BODY_CODE += ( -140120                   )

      NAIF_BODY_NAME += ( 'DIF_MRI'                 )
      NAIF_BODY_CODE += ( -140200                   )

   \begintext
      
   Deep Impact Impactor (DII) spacecraft and instruments IDs:

   \begindata

      NAIF_BODY_NAME += ( 'DEEP_IMPACT_IMPACTOR_SC' )
      NAIF_BODY_CODE += ( -70                       )

      NAIF_BODY_NAME += ( 'DII'                     )
      NAIF_BODY_CODE += ( -70                       )
      
      NAIF_BODY_NAME += ( 'DII_ITS'                 )
      NAIF_BODY_CODE += ( -70100                    )
      
   \begintext

   TEMPEL Comet ID:

   \begindata

      NAIF_BODY_NAME += ( 'TEMPEL'                  )
      NAIF_BODY_CODE += ( 1000093                   )
      
   \begintext


Deep Impact Frames
-------------------------------------------------------------------------------

   The following DIF and DII frames are defined in this kernel file:

           Name                  Relative to           Type       NAIF ID
      ======================  ===================  ============   =======

   Flyby Spacecraft Frame:
   -----------------------
      DIF_SPACECRAFT          rel.to J2000         CK             -140000

   Flyby Spacecraft Instrument Frames:
   -----------------------------------
      DIF_HRI_OPTICS          rel.to SPACECRAFT    FIXED          -140100
      DIF_HRI_VIS             rel.to SPACECRAFT    FIXED          -140110
      DIF_HRI_VIS_ADCS        rel.to HRI_VIS       FIXED          -140119
      DIF_HRI_IR              rel.to SPACECRAFT    FIXED          -140120
      DIF_MRI                 rel.to SPACECRAFT    FIXED          -140200
      DIF_MRI_ADCS            rel.to MRI           FIXED          -140209

   Flyby Spacecraft Star Tracker Frames:
   -------------------------------------
      DIF_ST1                 rel.to SPACECRAFT    FIXED          -140710
      DIF_ST2                 rel.to SPACECRAFT    FIXED          -140720

   Flyby Spacecraft Sun Sensor Frames:
   -----------------------------------
      DIF_CSS1                rel.to SPACECRAFT    FIXED          -140801
      DIF_CSS2                rel.to SPACECRAFT    FIXED          -140802
      DIF_CSS3                rel.to SPACECRAFT    FIXED          -140803
      DIF_CSS4                rel.to SPACECRAFT    FIXED          -140804
      DIF_CSS5                rel.to SPACECRAFT    FIXED          -140805
      DIF_CSS6                rel.to SPACECRAFT    FIXED          -140806
      DIF_CSS7                rel.to SPACECRAFT    FIXED          -140807
      DIF_CSS8                rel.to SPACECRAFT    FIXED          -140808
      DIF_CSS9                rel.to SPACECRAFT    FIXED          -140809 (*)
      DIF_CSS10               rel.to SPACECRAFT    FIXED          -140810
      DIF_CSS11               rel.to SPACECRAFT    FIXED          -140811
      DIF_CSS12               rel.to SPACECRAFT    FIXED          -140812
      DIF_CSS13               rel.to SPACECRAFT    FIXED          -140813

   Flyby Spacecraft Low Gain Antenna Frames:
   -----------------------------------------
      DIF_LGA-Y               rel.to SPACECRAFT    FIXED          -140910
      DIF_LGA+Y               rel.to SPACECRAFT    FIXED          -140920
      DIF_HGA_GIMBAL          rel.to SPACECRAFT    CK             -140930
      DIF_HGA                 rel.to HGA_GIMBAL    FIXED          -140931

   Impactor Spacecraft Frame:
   --------------------------
      DII_SPACECRAFT          rel.to J2000         CK             -70000

   Impactor Spacecraft Instrument Frame:
   -------------------------------------
      DII_ITS                 rel.to SPACECRAFT    FIXED          -70100
      DII_ITS_ADCS            rel.to ITS           FIXED          -70109

   Comet Frame:
   ------------
      TEMPEL_FIXED            rel.to J2000         PCK            1000093 


   (*) Although the frame for CSS9 is listed in the table, this FK file 
       does not contain a definition for it because CSS9 has been removed
       from the s/c per ECR from December 17, 2001


Deep Impact Frames Hierarchy
-------------------------------------------------------------------------------

   The diagram below shows Deep Impact frames hierarchy:


                      "J2000" INERTIAL
        +--------------------------------------------------------------+
        |                   |                     |                    .
        |<-pck              |                     |<-pck               .
        |                   |                     |                    .
        V                   |                     V                    .
    "TEMPEL_FIXED"          |                "IAU_EARTH"               .
    TEMPEL BFR(*)           |<-ck            EARTH BFR(*)              .
    -------------           |                ------------              .
                            |                                          .
                            |                                          .
                            |            DIF_HGA                       .
                            |            -------                       .
                            |               ^                          .
                            |               |                          .
                            |               |<-fixed                   .
                            |               |                          .
                            |        DIF_HGA_GIMBAL                    .
                            |        --------------                    .
                            |               ^                          .
                            |               |                          .
                            V               |<-ck                      .
                      "DIF_SPACECRAFT"      |                          .
        +------------------------------------------------------+       .
        |                   |               |      |      |    .       .
        |<-fixed            |<-fixed        |      |      |    .<-ck   .<-ck(#)
        |                   |               |      |      |    .       .
        |                   V               |      |      |    V       V
        |           "DIF_HRI_OPTICS"        |      |      |  "DII_SPACECRAFT"
        |           +---------------+       |      |      |  +--------------+
        |           |               |       |      |      |         |
        |    fixed->|        fixed->|       |      |      |  fixed->|
        |           |               |       |      |      |         |
        V            V              V       |      |      |         V
    "DIF_MRI"  "DIF_HRI_IR"  "DIF_HRI_VIS"  |      |      |     "DII_ITS"  
    ---------  ------------  -------------  |      |      |     ---------   
        |                           |       |      |      |         |
        |<-fixed             fixed->|       |      |      |  fixed->|
        |                           |       |      |      |         |
        V                           V       |      |      |         V
    "DIF_MRI_ADCS"      "DIF_HRI_VIS_ADCS"  |      |      |  "DII_ITS_ADCS"
    --------------      ------------------  |      |      |  --------------
                                            |      |      | 
                                     fixed->| fxd->|      |<-fixed 
                                            |      |      | 
                                            V      V      V
                                 DIF_LGA+Y/-Y DIF_ST(N) DIF_CSS(N)
                                 ------------ --------- ----------

   (*)  BFR -- body-fixed rotating frame

   (#)  DII orientation provided in CK is likely to be with respect to the 
        DIF prior to separation and with respect to an inertial frame after
        separation.


Flyby Spacecraft Bus Frame
-------------------------------------------------------------------------------
 
   The DIF spacecraft frame is defined by the s/c design as follows [from 4]:

      -  +X axis is aligned with the Launch Vehicle (LV);
 
      -  +Y axis is parallel to the solar array normal;

      -  +Z axis completes the right hand frame;

      -  the origin of the frame is at the LV separation plane, along the
         LV centerline.

   This diagram illustrates DIF s/c frame (based on Figure 3-2 from [4]):
                             
                              |\
                      ._______| \________________   
                      |       |  \      ||       \
                      |   HGA |  /|     ||        \
                      |       | / |     ||         \     Solar Array
                      |       |/ / \    ||          | (behind the s/c) 
                      | ._______/___\___.|          |
                      | |  _            ||          |
                      | .'   \          ||          |
                      .'  MRI |         ||          |
                    .' .'    /          ||          |
                  .' .'|    /   __      ||          |
       MRI        `.' ||___/  .'   \    ||          |
    Boresight   ,'    | |   .'      |   ||          |
              .'      | | .'   HRI  |   ||          |
             V        | .'         /    ||          |
                      .'   ,'|____/     ||          |
                    .'   ,'     ^ +Xsc  |.__________.
                  .'   .'______ | ______.
                 '   .'      |_ | _|
       HRI       .`.'       /__ | __\
    Boresight  .'               |
             .'         <-------x    . . . . . . . . . .  LV interface
            V       +Zsc        +Ysc                          plane
                            (into page) 


   The DIF s/c bus attitude is provided in CK, and, therefore, this frame is 
   defined as a CK-based frame.

   \begindata

      FRAME_DIF_SPACECRAFT     = -140000
      FRAME_-140000_NAME       = 'DIF_SPACECRAFT'
      FRAME_-140000_CLASS      = 3
      FRAME_-140000_CLASS_ID   = -140000
      FRAME_-140000_CENTER     = -140
      CK_-140000_SCLK          = -140
      CK_-140000_SPK           = -140

   \begintext


DIF Science Instrument Frames
-------------------------------------------------------------------------------

   This section contains frame definitions for DIF science instruments --
   High Resolution Imager (HRI) and Medium Resolution Imager (MRI.)


HRI Frames
----------

   HRI consists of a multi-spectral CCD and an infrared imaging
   spectrometer mounted on the same optical bench attached to the -Y
   side of the spacecraft bus.

   There are three "primary" frames defined for HRI: one for optical
   bench and one for each of the detectors. Nominally, all three frames
   are co-aligned. The diagram below illustrates the orientation of
   these frames with respect to the spacecraft:

                              |\
                      ._______| \________________   
                      |       |  \      ||       \
                      |   HGA |  /|     ||        \
                      |       | / |     ||         \     Solar Array
                      |       |/ / \    ||          | (behind the s/c) 
                      | ._______/___\___.|          |
                      | | __            ||          |
                      | .'              ||          |
                      .'    ^ +Xhri     ||          |
                    .' .'    \          ||          |
                  .' .'|      \ __                   
       MRI        `.' ||___/  .\  +Yhri (into the page)
    Boresight   ,'    | |   .' .x   |               
              .'          .' .'     |   ||          |
             V       +Zhri .'  HRI /    ||          |
                      .' v','|____/     ||          |
                    .'   ,'             ||__________.
                  .'   .'_______________.
                 '   .'      |_ ^ +Xsc
       HRI       .`.'       /__ | __\
    Boresight  .'               |
             .'            <----x    . . . . . . . . . .  LV interface
            V           +Zsc   +Ysc                          plane
                              (into page) 

   On this diagram +Xhri, +Yhri, +Zhri represent corresponding axes of 
   the DIF_HRI_OPTICS, DIF_HRI_VIS and DIF_HRI_IR frames.


HRI Optics Frame
----------------

   The axes of the DIF_HRI_OPTICS frame are defined by the instrument 
   design as follows:

      -  +Z axis is in the direction of HRI boresight; 
         
      -  +Y axis is parallel to the nominal HRI VIS CCD columns and HRI
         IR spatial resolution direction and nominally points in the
         same direction as the s/c +Y axis;
         
      -  +X completes the right hand frame;
      
      -  the origin of this frame is the HRI optics focal point.

   Nominally, HRI Optics frame is rotated with respect to the DIF spacecraft 
   frame by -45 degrees about spacecraft +Y axis: 

           hrio
          M    = | 0.0 |  * | -45.0 |  * | 0.0 |
           sc           Z            Y          X

   Note that this nominal alignment is valid during flight because
   instead of updating it after HRI to s/c calibration was performed,
   the star tracker alignments to s/c frame were adjusted to make sure
   that the nominal relationship between HRI and s/c frames is
   preserved.

   (The frame definition below contains the opposite of this rotation 
   because Euler angles specified in it define transformation from HRI 
   optics to s/c frame -- see [1].)

   \begindata

      FRAME_DIF_HRI_OPTICS      = -140100
      FRAME_-140100_NAME        = 'DIF_HRI_OPTICS'
      FRAME_-140100_CLASS       = 4
      FRAME_-140100_CLASS_ID    = -140100
      FRAME_-140100_CENTER      = -140
      TKFRAME_-140100_SPEC      = 'ANGLES'
      TKFRAME_-140100_RELATIVE  = 'DIF_SPACECRAFT'
      TKFRAME_-140100_ANGLES    = ( 0.0, 45.0, 0.0 )
      TKFRAME_-140100_AXES      = ( 1,    2,   3   )
      TKFRAME_-140100_UNITS     = 'DEGREES'

   \begintext


HRI VIS and IR Frames
---------------------

   The axes of the DIF_HRI_VIS frame are defined along the lines of the
   standard image frame convention:

      -  +Z axis is along the instrument boresight;
         
      -  +X axis is along the instrument CCD lines and point toward the
         right edge of the image;
         
      -  +Y axis is along the instrument CCD columns; it completes the
         right hand frame and points from toward the bottom of the
         image;
      
      -  the origin of this frame is located at the geometric center of the 
         instrument CCD.

   The axes of the DIF_HRI_IR frame are defined as follows:

      -  +Z axis is along the instrument boresight;
         
      -  +Y axis is along the spatial resolution direction of the 
         instrument's CCD;

      -  +X axis is along the spectral resolution direction of the 
         instrument CCD; it completes the right hand frame;
         
      -  the origin of this frame is located at the geometric center of
         the instrument CCD.
  
   Nominally both frames, HRI VIS and HRI IR, are co-aligned with the 
   HRI optics frame.

   While the zero nominal alignment with respect to the HRI optics
   frame is valid for DIF_HRI_VIS frame during flight (see note in the
   ``HRI Optics Section'' above), it is not valid for the DIF_HRI_IR
   frame. The actual alignment the HRI IR frame was derived from the
   flight data by the DI science team and given in [12]. According this
   solution the HRI IR boresight -- +Z axis of the DIF_HRI_IR frame
   corresponding to the center of stored full-frame IR image -- is
   tilted with respect to the HRI VIS boresight -- +Z axis of the
   DIF_HRI_VIS frame -- by 36 microradians toward -X and by 7
   microradians toward +Y. Since these two angles have very small
   magnitude they are treated as rotation angles and incorporated ``as
   is'' into the DIF_HRI_IR frame definition below. The third angle is
   set to zero because the same analysis determined that the rotation
   of the HRI IR slit about the boresight with respect to the MRI/HRI
   is very small, sub-pixel level and for this reason could be
   neglected.

   (The frame definitions below contain the opposite of this rotation 
   because Euler angles specified in it define transformation from HRI 
   IR/VIS to optics frame -- see [1].)

   \begindata

      FRAME_DIF_HRI_VIS         = -140110
      FRAME_-140110_NAME        = 'DIF_HRI_VIS'
      FRAME_-140110_CLASS       = 4
      FRAME_-140110_CLASS_ID    = -140110
      FRAME_-140110_CENTER      = -140
      TKFRAME_-140110_SPEC      = 'ANGLES'
      TKFRAME_-140110_RELATIVE  = 'DIF_HRI_OPTICS'
      TKFRAME_-140110_ANGLES    = ( 0.0, 0.0, 0.0 )
      TKFRAME_-140110_AXES      = ( 1,   2,   3  )
      TKFRAME_-140110_UNITS     = 'DEGREES'

      FRAME_DIF_HRI_IR          = -140120
      FRAME_-140120_NAME        = 'DIF_HRI_IR'
      FRAME_-140120_CLASS       = 4
      FRAME_-140120_CLASS_ID    = -140120
      FRAME_-140120_CENTER      = -140
      TKFRAME_-140120_SPEC      = 'ANGLES'
      TKFRAME_-140120_RELATIVE  = 'DIF_HRI_OPTICS'
      TKFRAME_-140120_ANGLES    = ( 0.0, 0.000036, 0.000007 )
      TKFRAME_-140120_AXES      = ( 3,   2,        1        )
      TKFRAME_-140120_UNITS     = 'RADIANS'

   \begintext


HRI VIS ADSC Frame
------------------

   The "primary" HRI VIS frame (DIF_HRI_VIS) defined above follows the
   OPNAV convention. The ADCS defined another frame for this camera
   using a different convention. The axes of the ADCS frame, in this FK
   named DIF_HRI_VIS_ADCS, are related to the axes of the "primary" HRI
   VIS frame (DIF_HRI_VIS) as follows (see [8]):
   
      -  +X axis of DIF_HRI_VIS_ADCS is the same as +Z axis of
         DIF_HRI_VIS (along boresight);

      -  +Y axis of DIF_HRI_VIS_ADCS is the same as +Y axis of
         DIF_HRI_VIS (along CCD columns);

      -  +Z axis of DIF_HRI_VIS_ADCS is the same as -X axis of
         DIF_HRI_VIS (along CCD lines);

   A single rotation of -90 degrees about +Y is needed to align the 
   DIF_HRI_VIS with DIF_HRI_VIS_ADCS frame.

   (The frame definitions below contain the opposite of this rotation 
   because Euler angles specified in it define transformation from ADCS 
   to "primary" frame -- see [1].)

   \begindata

      FRAME_DIF_HRI_VIS_ADCS    = -140119
      FRAME_-140119_NAME        = 'DIF_HRI_VIS_ADCS'
      FRAME_-140119_CLASS       = 4
      FRAME_-140119_CLASS_ID    = -140119
      FRAME_-140119_CENTER      = -140
      TKFRAME_-140119_SPEC      = 'ANGLES'
      TKFRAME_-140119_RELATIVE  = 'DIF_HRI_VIS'
      TKFRAME_-140119_ANGLES    = ( 0.0, 90.0, 0.0 )
      TKFRAME_-140119_AXES      = ( 1,    2,   3   )
      TKFRAME_-140119_UNITS     = 'DEGREES'

   \begintext


MRI Frame
---------

   MRI consists of a single multi-spectral CCD mounted on an optical
   bench attached to the -Y side of the spacecraft bus.

   Only one frame is defined for MRI. This frame orientation with
   respect to the spacecraft is shown on the diagram below:

                              |\
                      ._______| \________________   
                      |       |  \      ||       \
                          HGA |  /|     ||        \
                      ^       | / |     ||         \     Solar Array
                 +Xmri \      |/ / \    ||          |  (behind the s/c) 
                      | \  _____/___\___.|          |
                      |  \__                        |
                      | .'x +Ymri (into the page )  |
                      . .'                          |
                +Zmri'.'.'   /          ||          |
                  .' V.'    /   __      ||          |             
       MRI        ` .'||___/  .'   \    ||          |
    Boresight   ,'    | |   .'      |   ||          |            
              .'          .'        |   ||          |
             V          .'     HRI /    ||          |
                      .'   ,'|____/     ||          |
                    .'   ,'             ||__________.
                  .'   .'_______________.
                 '   .'      |_ ^ +Xsc
       HRI       .`.'       /__ | __\
    Boresight  .'               |
             .'            <----x    . . . . . . . . . .  LV interface
            V           +Zsc   +Ysc                          plane
                              (into page) 


   The axes of the DIF_MRI frame are defined by the instrument design 
   as follows:

      -  +Z axis is in the direction of MRI boresight; 
         
      -  +Y axis is parallel to the nominal MRI CCD columns and points
         toward the bottom of the image; it nominally points in the
         same direction as the s/c +Y axis;
         
      -  +X completes the right hand frame; it is parallel to MRI CCD
         lines and points toward the right side of the image;
      
      -  the origin of this frame is the MRI optics focal point.

   Nominally, MRI frame is rotated with respect to the DIF spacecraft 
   frame by -45 degrees about spacecraft +Y axis: 

           mri
          M    = | 0.0 |  * | -45.0 |  * | 0.0 |
           sc           Z            Y          X

   Initial estimate of the actual MRI alignment relative to the s/c was
   performed by DI OPNAV using images taken on 1-14-05. The following
   quaternion -- non-SPICE style, rotating from s/c to MRI -- for
   this misalignment was provided by Nick Mastrodemos, DI OPNAV ([10]):

          qb2mri   = [-0.00040450121923 -0.38272258737675   0.00104114520178   
                       0.92386263779150];

   Rotation angles corresponding to this quaternion are:

           mri
          M    = |0.130801483771|  * |-45.004882929701|  * |0.004013880181|
           sc                    Z                     Y                   X
   
   Final estimate of the actual MRI alignment relative to the s/c was
   calculated by DI AACS using images taken on 2005-108 and 2005-113.
   The following quaternion -- non-SPICE style, rotating from s/c to
   MRI -- for this misalignment was provided in [11]:

          qb2mri   = [-0.00039317053464, -0.38273873987874, 0.00102800577131, 
                       0.92385596584260];

   Rotation angles corresponding to this quaternion are:


           mri
          M    = |0.129539306414|  * |-45.006884881185|  * |0.004898709285|
           sc                    Z                     Y                   X
   
   (The frame definition below contains the opposite of this rotation 
   because Euler angles specified in it define transformation from MRI 
   to s/c frame -- see [1].)

   \begindata

      FRAME_DIF_MRI             = -140200
      FRAME_-140200_NAME        = 'DIF_MRI'
      FRAME_-140200_CLASS       = 4
      FRAME_-140200_CLASS_ID    = -140200
      FRAME_-140200_CENTER      = -140
      TKFRAME_-140200_SPEC      = 'ANGLES'
      TKFRAME_-140200_RELATIVE  = 'DIF_SPACECRAFT'
      TKFRAME_-140200_ANGLES    = ( -0.004898709285, 
                                    45.006884881185, 
                                    -0.129539306414 )
      TKFRAME_-140200_AXES      = ( 1,    2,   3   )
      TKFRAME_-140200_UNITS     = 'DEGREES'

   \begintext


MRI ADSC Frame
--------------

   The "primary" MRI frame (DIF_MRI) defined above follows the OPNAV
   convention. The ADCS defined another frame for this camera using a
   different convention. The axes of the ADCS frame, in this FK named
   DIF_MRI_ADCS, are related to the axes of the "primary" MRI frame
   (DIF_MRI) as follows (see [8]):
   
      -  +X axis of DIF_MRI_ADCS is the same as +Z axis of
         DIF_MRI (along boresight);

      -  +Y axis of DIF_MRI_ADCS is the same as +Y axis of
         DIF_MRI (along CCD columns);

      -  +Z axis of DIF_MRI_ADCS is the same as -X axis of
         DIF_MRI (along CCD lines);

   A single rotation of -90 degrees about +Y is needed to align the 
   DIF_MRI with DIF_MRI_ADCS frame.

   (The frame definitions below contain the opposite of this rotation 
   because Euler angles specified in it define transformation from ADCS 
   to "primary" frame -- see [1].)

   \begindata

      FRAME_DIF_MRI_ADCS        = -140209
      FRAME_-140209_NAME        = 'DIF_MRI_ADCS'
      FRAME_-140209_CLASS       = 4
      FRAME_-140209_CLASS_ID    = -140209
      FRAME_-140209_CENTER      = -140
      TKFRAME_-140209_SPEC      = 'ANGLES'
      TKFRAME_-140209_RELATIVE  = 'DIF_MRI'
      TKFRAME_-140209_ANGLES    = ( 0.0, 90.0, 0.0 )
      TKFRAME_-140209_AXES      = ( 1,    2,   3   )
      TKFRAME_-140209_UNITS     = 'DEGREES'

   \begintext


DIF Star Tracker and Sun Sensor Frames
-------------------------------------------------------------------------------

   This section contains frame definitions for DIF Star Trackers (ST1 and
   ST2) and Coarse Sun Sensors (CSS1 - CSS13).


Star Tracker Frames
-------------------

   The following description of the Star Tracker 1 and 2 (ST1 and ST2) 
   mounting alignment is provided in [5]:

      -- Tracker boresights are tilted 18 degrees above the instrument 
         platform (toward spacecraft -Y);

      -- Tracker 1 (ST1) is clocked 20 degrees of the HRI anti-boresight,
         toward the SIM;

      -- Tracker 2 (ST2) is clocked zero degrees of the HRI anti-boresight;

   This diagram illustrates ST mounting alignment and boresight directions:


                              |\
                      ._______| \________________   
                      |       |  \      ||    .> \
                      |   HGA |  /|     ||  .' ST2 boresight
                      |       | / |  ST2  .'       \ 
                 +Ymri|       |/ / \   .`.          |  
                      |    _____/___\.' .'        .> ST1 boresight
                      |   __       .' .'||     .-' 
                      | .'o  \      ''  ||.o'`.     |
                      .'  MRI |        .o' .o'      |
                    .'  .'   /         '.o'ST1      |
                  .'  .'    /   __      ||          |             
       MRI        ` .'||___/  .'   \    ||          |
    Boresight   ,'    | |   .'      |   ||          |            
              .'          .'        |   ||          |
             V          .'     HRI /    ||          |
                      .'   ,'|____/     ||          |    Solar Array
                    .'   ,'             ||__________.  (behind the s/c) 
                  .'   .'_______________.
                 '   .'      |_ ^ +Xsc
       HRI       .`.'       /__ | __\
    Boresight  .'               |
             .'            <----x    . . . . . . . . . .  LV interface
            V           +Zsc   +Ysc                          plane
                              (into page) 

      
   The frame for each of the two trackers is defined as follows:

      -  +Z axis is along the tracker's boresight;
         
      -  +X axis is parallel to the XZ plane of the DIF s/c frame;
         
      -  +Y axis completes the right hand frame and points approximately 
         in the direction of the s/c -Y;
      
      -  the origin of the frame is located at the geometric center of the 
         tracker.

   The following ST1 and ST2 boresight directions with respect to the 
   DIF s/c frame and spacecraft to tracker rotation matrices are provided 
   in [5]:

      Star Tracker boresight orientation:

         ST1:  [+0.401934, -0.309017, -0.861950]
  
         ST2:  [+0.672499, -0.309017, -0.672499]

      Spacecraft to tracker rotation matrix:

         ST1:    0.90630779          0.00000000          0.42261826
                -0.13059623         -0.95105652          0.28006451
                 0.40193385         -0.30901699         -0.86194993

         ST2:    0.70710678          0.00000000          0.70710678
                -0.21850801         -0.95105652          0.21850801
                 0.67249851         -0.30901699         -0.67249851

   Based on these directions/matrices and description of the mounting 
   the following rotations are need to transform the s/c frame into 
   the tracker frames:

      ST1:  first, by 155 degrees about +Y; second, by 18 degrees
            about +X; third, by 180 degrees about +Z; 
       
      ST2:  first, by 135 degrees about +Y; second, by 18 degrees 
            about +X; third, by 180 degrees about +Z;

   (The frame definitions below contains the opposite of this rotation 
   because Euler angles specified in it define transformation from tracker 
   to s/c frame -- see [1].)

   \begindata

      FRAME_DIF_ST1             = -140710
      FRAME_-140710_NAME        = 'DIF_ST1'
      FRAME_-140710_CLASS       = 4
      FRAME_-140710_CLASS_ID    = -140710
      FRAME_-140710_CENTER      = -140
      TKFRAME_-140710_SPEC      = 'ANGLES'
      TKFRAME_-140710_RELATIVE  = 'DIF_SPACECRAFT'
      TKFRAME_-140710_ANGLES    = ( -155.0, -18.0, -180.0 )
      TKFRAME_-140710_AXES      = (    2,     1,      3   )
      TKFRAME_-140710_UNITS     = 'DEGREES'

      FRAME_DIF_ST2             = -140720
      FRAME_-140720_NAME        = 'DIF_ST2'
      FRAME_-140720_CLASS       = 4
      FRAME_-140720_CLASS_ID    = -140720
      FRAME_-140720_CENTER      = -140
      TKFRAME_-140720_SPEC      = 'ANGLES'
      TKFRAME_-140720_RELATIVE  = 'DIF_SPACECRAFT'
      TKFRAME_-140720_ANGLES    = ( -135.0, -18.0, -180.0 )
      TKFRAME_-140720_AXES      = (    2,     1,      3   )
      TKFRAME_-140720_UNITS     = 'DEGREES'

   \begintext


Coarse Sun Sensor Frames
------------------------

   The frame for each of the 13 Coarse Sun Sensors (CSS's) mounted on the
   DIF is defined as follows:

      -  +Z axis is along the sensor's boresight;
         
      -  +Y axis is parallel to the XY plane of the DIF s/c frame;
         
      -  +X axis is completes the right hand frame;
      
      -  the origin of the frame is located at the geometric center of the 
         sensor.

   The following CSS boresight directions with respect to the DIF s/c frame
   are provided in the table 2 of [6]:

      Sensor   FOV, deg      X         Y         Z
      -------  --------- --------- --------- ---------
         1        53       0.6691    0.7431    0
         2        78      -0.8160    0.5714    0.0872
         3        78      -0.1228    0.6964   -0.7071
         4        78      -0.1228    0.6964    0.7071
         5        78       0.2743    0.9262    0.2588
         6        78       0.5649    0.0996   -0.8192
         7        78       0.5649    0.0966    0.8192
         8        78       0.0872   -0.9962   -0.0003
         9         -        -         -         -     (*)
         10       78       0.0517   -0.4492    0.8919
         11       78      -0.9513   -0.1677    0.2588
         12       78      -0.9513   -0.1677   -0.2588
         13       78      -0.9962    0.0872    0

   (*) CSS9 has been removed per ECR from December 17, 2001

   These directions can presented by the following azimuth/elevation pairs
   with azimuth measured CCW from +X in the s/c XY plane and elevation 
   measured positive from the XY plane toward +Z axis (azimuth and elevation 
   are defined consistently with the diagram on page 1 of [6]):

      Sensor       AZ        EL
      -------  --------- ---------
         1        48.0       0.0
         2       145.0       5.0
         3       100.0     -45.0
         4       100.0      45.0
         5        73.5      15.0
         6        10.0     -55.0
         7         9.7      55.0
         8       -85.0      -0.0
         9          -         -
         10      -83.4      63.1
         11     -170.0      15.0
         12     -170.0     -15.0
         13      175.0       0.0

   A matrix rotating vectors from the s/c frame to the sensor frame can 
   be constructed from the azimuth and elevation pair as follows:
  
         M = |90.0-EL|  * |AZ|
                      Y       Z

   (The frame definitions below contain the opposite of this rotation 
   because Euler angles specified in it define transformation from 
   sensor to s/c frame -- see [1].)

   \begindata

      FRAME_DIF_CSS1             = -140801
      FRAME_-140801_NAME         = 'DIF_CSS1'
      FRAME_-140801_CLASS        = 4
      FRAME_-140801_CLASS_ID     = -140801
      FRAME_-140801_CENTER       = -140
      TKFRAME_-140801_SPEC       = 'ANGLES'
      TKFRAME_-140801_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140801_ANGLES     = (    -48.000,    -90.000,      0.000 )
      TKFRAME_-140801_AXES       = (      3,          2,          1     )
      TKFRAME_-140801_UNITS      = 'DEGREES'

      FRAME_DIF_CSS2             = -140802
      FRAME_-140802_NAME         = 'DIF_CSS2'
      FRAME_-140802_CLASS        = 4
      FRAME_-140802_CLASS_ID     = -140802
      FRAME_-140802_CENTER       = -140
      TKFRAME_-140802_SPEC       = 'ANGLES'
      TKFRAME_-140802_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140802_ANGLES     = (   -145.000,    -85.000,      0.000 )
      TKFRAME_-140802_AXES       = (      3,          2,          1     )
      TKFRAME_-140802_UNITS      = 'DEGREES'

      FRAME_DIF_CSS3             = -140803
      FRAME_-140803_NAME         = 'DIF_CSS3'
      FRAME_-140803_CLASS        = 4
      FRAME_-140803_CLASS_ID     = -140803
      FRAME_-140803_CENTER       = -140
      TKFRAME_-140803_SPEC       = 'ANGLES'
      TKFRAME_-140803_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140803_ANGLES     = (   -100.000,   -135.000,      0.000 )
      TKFRAME_-140803_AXES       = (      3,          2,          1     )
      TKFRAME_-140803_UNITS      = 'DEGREES'

      FRAME_DIF_CSS4             = -140804
      FRAME_-140804_NAME         = 'DIF_CSS4'
      FRAME_-140804_CLASS        = 4
      FRAME_-140804_CLASS_ID     = -140804
      FRAME_-140804_CENTER       = -140
      TKFRAME_-140804_SPEC       = 'ANGLES'
      TKFRAME_-140804_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140804_ANGLES     = (   -100.000,    -45.000,      0.000 )
      TKFRAME_-140804_AXES       = (      3,          2,          1     )
      TKFRAME_-140804_UNITS      = 'DEGREES'

      FRAME_DIF_CSS5             = -140805
      FRAME_-140805_NAME         = 'DIF_CSS5'
      FRAME_-140805_CLASS        = 4
      FRAME_-140805_CLASS_ID     = -140805
      FRAME_-140805_CENTER       = -140
      TKFRAME_-140805_SPEC       = 'ANGLES'
      TKFRAME_-140805_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140805_ANGLES     = (    -73.500,    -75.000,      0.000 )
      TKFRAME_-140805_AXES       = (      3,          2,          1     )
      TKFRAME_-140805_UNITS      = 'DEGREES'

      FRAME_DIF_CSS6             = -140806
      FRAME_-140806_NAME         = 'DIF_CSS6'
      FRAME_-140806_CLASS        = 4
      FRAME_-140806_CLASS_ID     = -140806
      FRAME_-140806_CENTER       = -140
      TKFRAME_-140806_SPEC       = 'ANGLES'
      TKFRAME_-140806_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140806_ANGLES     = (    -10.000,   -145.000,      0.000 )
      TKFRAME_-140806_AXES       = (      3,          2,          1     )
      TKFRAME_-140806_UNITS      = 'DEGREES'

      FRAME_DIF_CSS7             = -140807
      FRAME_-140807_NAME         = 'DIF_CSS7'
      FRAME_-140807_CLASS        = 4
      FRAME_-140807_CLASS_ID     = -140807
      FRAME_-140807_CENTER       = -140
      TKFRAME_-140807_SPEC       = 'ANGLES'
      TKFRAME_-140807_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140807_ANGLES     = (     -9.700,    -35.000,      0.000 )
      TKFRAME_-140807_AXES       = (      3,          2,          1     )
      TKFRAME_-140807_UNITS      = 'DEGREES'

      FRAME_DIF_CSS8             = -140808
      FRAME_-140808_NAME         = 'DIF_CSS8'
      FRAME_-140808_CLASS        = 4
      FRAME_-140808_CLASS_ID     = -140808
      FRAME_-140808_CENTER       = -140
      TKFRAME_-140808_SPEC       = 'ANGLES'
      TKFRAME_-140808_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140808_ANGLES     = (     85.000,    -90.000,      0.000 )
      TKFRAME_-140808_AXES       = (      3,          2,          1     )
      TKFRAME_-140808_UNITS      = 'DEGREES'

      FRAME_DIF_CSS10            = -140810
      FRAME_-140810_NAME         = 'DIF_CSS10'
      FRAME_-140810_CLASS        = 4
      FRAME_-140810_CLASS_ID     = -140810
      FRAME_-140810_CENTER       = -140
      TKFRAME_-140810_SPEC       = 'ANGLES'
      TKFRAME_-140810_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140810_ANGLES     = (     83.400,    -26.900,      0.000 )
      TKFRAME_-140810_AXES       = (      3,          2,          1     )
      TKFRAME_-140810_UNITS      = 'DEGREES'

      FRAME_DIF_CSS11            = -140811
      FRAME_-140811_NAME         = 'DIF_CSS11'
      FRAME_-140811_CLASS        = 4
      FRAME_-140811_CLASS_ID     = -140811
      FRAME_-140811_CENTER       = -140
      TKFRAME_-140811_SPEC       = 'ANGLES'
      TKFRAME_-140811_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140811_ANGLES     = (    170.000,    -75.000,      0.000 )
      TKFRAME_-140811_AXES       = (      3,          2,          1     )
      TKFRAME_-140811_UNITS      = 'DEGREES'

      FRAME_DIF_CSS12            = -140812
      FRAME_-140812_NAME         = 'DIF_CSS12'
      FRAME_-140812_CLASS        = 4
      FRAME_-140812_CLASS_ID     = -140812
      FRAME_-140812_CENTER       = -140
      TKFRAME_-140812_SPEC       = 'ANGLES'
      TKFRAME_-140812_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140812_ANGLES     = (    170.000,   -105.000,      0.000 )
      TKFRAME_-140812_AXES       = (      3,          2,          1     )
      TKFRAME_-140812_UNITS      = 'DEGREES'

      FRAME_DIF_CSS13            = -140813
      FRAME_-140813_NAME         = 'DIF_CSS13'
      FRAME_-140813_CLASS        = 4
      FRAME_-140813_CLASS_ID     = -140813
      FRAME_-140813_CENTER       = -140
      TKFRAME_-140813_SPEC       = 'ANGLES'
      TKFRAME_-140813_RELATIVE   = 'DIF_SPACECRAFT'
      TKFRAME_-140813_ANGLES     = (   -175.000,    -90.000,      0.000 )
      TKFRAME_-140813_AXES       = (      3,          2,          1     )
      TKFRAME_-140813_UNITS      = 'DEGREES'

   \begintext


DIF Antenna Frames
-------------------------------------------------------------------------------

   This section contains frame definitions for DIF Low Gain (LGA+Y and
   LGA-Y) and High Gain (HGA) Antenna Frames.


DIF Low Gain Antenna Frames
---------------------------

   This section contains frame definitions for DIF Low Gain Antennas, 
   LGA-Y and LGA+Y.

   This diagram illustrates LGA mounting alignment and boresight 
   directions (from [7]):


                            |\  HGA
                      ..    | \
                      ||    |  \ 
                      ||    |  /|
            +Xlga+y ^ ||    | / |             ^ +Xlga-y
                    | ||    |/ / \            |
                    | ||._____/___\_____.     |
       +Zlga+y <----x#|||               |====#o----> +Zlga-y
              +Ylga+y |||               |-. +Ylga-y
               (into  |||               | | (out of
                page) |||               | |  page)
                      |||               | |
                      |||               | ||
                      |||               | ||
                      |||               |_.| HRI
                      |||               |  |
                      |||               |  |
               Solar  |||               |__|
               Array    ._______________.    Science Deck
                             |_ ^ +Xsc
                            /__ | __\
                      +Ysc      |     
                           <----o     . . . . . . . . .  LV interface
                               +Zsc                           plane
                              (out of
                               page) 

      
   The frame for each of the antennas is defined as follows:

      -  +Z axis is along the antenna's boresight;
         
      -  +X axis is along the antenna pattern reference axis;
         
      -  +Y axis completes the right hand frame;
      
      -  the origin of the frame is located at the geometric center of the 
         antenna patch.

   Assuming that the antenna pattern reference direction is along the
   s/c +X axis, a single rotation by +90 degrees about +X axis is needed to 
   align the s/c frame with the LGA-Y frame and a single rotation by -90 
   degrees about +X axis is needed to align the s/c frame with the LGA+Y 
   frame.

   (The frame definitions below contains the opposite of this rotation 
   because Euler angles specified in it define transformation from tracker 
   to s/c frame -- see [1].)

   \begindata

      FRAME_DIF_LGA-Y           = -140910
      FRAME_-140910_NAME        = 'DIF_LGA-Y'
      FRAME_-140910_CLASS       = 4
      FRAME_-140910_CLASS_ID    = -140910
      FRAME_-140910_CENTER      = -140
      TKFRAME_-140910_SPEC      = 'ANGLES'
      TKFRAME_-140910_RELATIVE  = 'DIF_SPACECRAFT'
      TKFRAME_-140910_ANGLES    = (    0.0,   0.0,  -90.0 )
      TKFRAME_-140910_AXES      = (    3,     2,      1   )
      TKFRAME_-140910_UNITS     = 'DEGREES'

      FRAME_DIF_LGA+Y           = -140920
      FRAME_-140920_NAME        = 'DIF_LGA+Y'
      FRAME_-140920_CLASS       = 4
      FRAME_-140920_CLASS_ID    = -140920
      FRAME_-140920_CENTER      = -140
      TKFRAME_-140920_SPEC      = 'ANGLES'
      TKFRAME_-140920_RELATIVE  = 'DIF_SPACECRAFT'
      TKFRAME_-140920_ANGLES    = (    0.0,   0.0,   90.0 )
      TKFRAME_-140920_AXES      = (    3,     2,      1   )
      TKFRAME_-140920_UNITS     = 'DEGREES'

   \begintext


DIF High Gain Antenna Frames
----------------------------

   This section contains frame definitions for DIF High Gain Antenna 
   (HGA) frames. Two frames are defined for HGA -- HGA gimbal frame, 
   DIF_HGA_GIMBAL, and HGA frame, DIF_HGA.

   The HGA gimbal frame, DIF_HGA_GIMBAL, is defined as follows:

      -  +X axis is along the antenna's boresight;
         
      -  +Y axis is along the HGA gimbal ``inboard'' rotation axis;
         
      -  +Z axis is along the HGA gimbal ``outboard'' rotation axis;
      
      -  the origin of the frame is located at the intersection of 
         the antenna ``outboard'' gimbal axis and antenna boresight.

   The HGA frame, DIF_HGA, is defined as follows:

      -  +Z axis is along the antenna's boresight;
         
      -  +X axis is along the HGA gimbal ``outboard'' rotation axis and
         points in the direction opposite to the HGA gimbal +Z axis;
         
      -  +Y axis completes the right hand frame;
      
      -  the origin of the frame is located at the intersection of 
         the antenna boresight and outer rim plane of the antenna dish.

   This diagram illustrates the two HGA frame (HGA is shown in zero 
   gimbal position):

                                ^ +Zhga
                                |                     
                         +Yhga  |                      +Xhga is into
                           <----x_____                    the page.
                           `.   ^ +Xhgag
                      ..     `. | .'
                      ||   <----o'                    +Zhgag is out of 
                      || +Yhgag |                         the page.
                      ||        |             
                      ||       / \            
                      ||._____/___\_____.     
                     #|||               |====#
                      |||               |-.
                      |||               | |
                      |||               | |
                      |||               | |
                      |||               | ||
                      |||               | ||
                      |||               |_.| HRI
                      |||               |  |
                      |||               |  |
               Solar  |||               |__|
               Array    ._______________.    Science Deck
                             |_ ^ +Xsc
                            /__ | __\
                      +Ysc      |     
                           <----o     . . . . . . . . .  LV interface
                               +Zsc                           plane
                           (out of page) 



   The HGA gimbal frame is defined as a CK-based frame because its
   orientation with respect to the s/c, determined by the rotations in
   ``inboard'' (first, about Y) and ``outboard'' (second, about Z), is
   stored in CK files. In zero gimbal position this frame is
   co-aligned with the s/c frame.

   The HGA frame is defined as a fixed offset frame with respect to the
   HGA gimbal frame. The main purpose of introducing this frame is to
   have a frame in which the antenna boresight is co-aligned with the
   +Z axis of the frame. A single rotation by +90 degrees about +Y axis
   is needed to lineup the HGA gimbal frame with the HGA frame.
   
   \begindata

      FRAME_DIF_HGA_GIMBAL     = -140930
      FRAME_-140930_NAME       = 'DIF_HGA_GIMBAL'
      FRAME_-140930_CLASS      = 3
      FRAME_-140930_CLASS_ID   = -140930
      FRAME_-140930_CENTER     = -140
      CK_-140930_SCLK          = -140
      CK_-140930_SPK           = -140

      FRAME_DIF_HGA            = -140931
      FRAME_-140931_NAME       = 'DIF_HGA'
      FRAME_-140931_CLASS      = 4
      FRAME_-140931_CLASS_ID   = -140931
      FRAME_-140931_CENTER     = -140
      TKFRAME_-140931_SPEC     = 'ANGLES'
      TKFRAME_-140931_RELATIVE = 'DIF_HGA_GIMBAL'
      TKFRAME_-140931_ANGLES   = (    0.0,   0.0,  -90.0 )
      TKFRAME_-140931_AXES     = (    3,     1,      2   )
      TKFRAME_-140931_UNITS    = 'DEGREES'

   \begintext


Impactor Spacecraft Bus Frame
-------------------------------------------------------------------------------
 
   The DII spacecraft frame is defined by the s/c design as follows:

      -  +X axis points in the direction opposite to the ITS boresight;
 
      -  +Y axis points in the direction of the combined thrust of the
         two side-mounted attitude control thrusters;

      -  +Z axis completes the right hand frame;

      -  the origin of the frame is at launch vehicle interface.

   This diagram illustrates DII s/c frame:

                                 ._____.    ^ +Zsc
                   Star        .'       `..'
                 Tracker   .`./         .'\
                           ___`.      .'   .
                             |      o'     |
                             .  +Xsc `.    @ Thruster
                              \        `. /
                               `.       .`.
                                 ` ---@'   `> +Ysc
                                 Thruster            +Xsc is out of 
                                                         the page.

   The DII s/c bus attitude is provided in CK, and, therefore, this
   frame is defined as a CK-based frame.

   \begindata

      FRAME_DII_SPACECRAFT     = -70000
      FRAME_-70000_NAME        = 'DII_SPACECRAFT'
      FRAME_-70000_CLASS       = 3
      FRAME_-70000_CLASS_ID    = -70000
      FRAME_-70000_CENTER      = -70
      CK_-70000_SCLK           = -70
      CK_-70000_SPK            = -70

   \begintext


DII Science Instrument Frames
-------------------------------------------------------------------------------

   This section contains frame definitions for DII science instrument
   -- Impactor Targeting Sensor (ITS)


ITS Frame
---------

   ITS consists of a single multi-spectral CCD mounted on an optical
   bench attached to the impactor spacecraft bus.

   Only one frame is defined for ITS. This frame orientation with
   respect to the DII spacecraft is shown on the diagram below:

                                 ._____.    ^ +Zsc
                   Star        +Zits    `..'
                 Tracker   .`./     x--------> +Xits
                           ___`.    | .'   .
                             |      o'     |
                             . +Xsc |`.    @ Thruster
                              \     |  `. /
                              +Yits V   .`.
                                 ` ---@'   `> +Ysc
                                       Thruster 
                                                   +Xsc is out of the page.
                                                   +Zits is into the page.

   The axes of the DIF_ITS frame are defined by the instrument design 
   as follows:

      -  +Z axis is in the direction of ITS boresight; it nominally 
         points in the direction of the DII -X axis;
         
      -  +X axis is parallel to ITS CCD lines and points toward the 
         right side of the image; it is nominally 30 degrees off the 
         s/c +Y axis toward the s/c +Z axis;

      -  +Y axis is parallel to the nominal ITS CCD columns and points
         toward the bottom of the image; it completes the right hand
         frame;
      
      -  the origin of this frame is the ITS optics focal point.

   Nominally, ITS frame is rotated with respect to the DII spacecraft
   frame by -90 degrees about spacecraft +Y axis and then by +60
   degrees about +Z axis:

           its
          M    = | +60.0 |  * | -90.0 |  * | 0.0 |
           sc             Z            Y          X

   (The frame definition below contains the opposite of this rotation 
   because Euler angles specified in it define transformation from ITS 
   to s/c frame -- see [1].)

   \begindata

      FRAME_DII_ITS             = -70100
      FRAME_-70100_NAME         = 'DII_ITS'
      FRAME_-70100_CLASS        = 4
      FRAME_-70100_CLASS_ID     = -70100
      FRAME_-70100_CENTER       = -70
      TKFRAME_-70100_SPEC       = 'ANGLES'
      TKFRAME_-70100_RELATIVE   = 'DII_SPACECRAFT'
      TKFRAME_-70100_ANGLES     = ( 0.0, 90.0, -60.0 )
      TKFRAME_-70100_AXES       = ( 1,    2,     3   )
      TKFRAME_-70100_UNITS      = 'DEGREES'

   \begintext


ITS ADSC Frame
--------------

   The "primary" ITS frame (DII_ITS) defined above follows the OPNAV
   convention. The ADCS defined another frame for this camera using a
   different convention. The axes of the ADCS frame, in this FK named
   DII_ITS_ADCS, are related to the axes of the "primary" ITS frame
   (DII_ITS) and Impactor spacecraft bus frame (DII_SPACECRAFT) as
   follows (see [8]):
   
      -  +X axis of DII_ITS_ADCS is the same as +Z axis of
         DII_ITS (along boresight);

      -  +Y axis of DII_ITS_ADCS is the same as +Y axis of
         DII_SPACECRAFT;

      -  +Z axis of DII_ITS_ADCS is the same as -Z axis of
         DII_SPACECRAFT;

   Two rotations are needed to align the DII_ITS with DII_ITS_ADCS
   frame: first by -90 degrees about +Y, then by -60 degrees about +X.

   (The frame definitions below contain the opposite of this rotation
   because Euler angles specified in it define transformation from ADCS
   to "primary" frame -- see [1].)

   \begindata

      FRAME_DII_ITS_ADCS        = -70109
      FRAME_-70109_NAME         = 'DII_ITS_ADCS'
      FRAME_-70109_CLASS        = 4
      FRAME_-70109_CLASS_ID     = -70109
      FRAME_-70109_CENTER       = -70
      TKFRAME_-70109_SPEC       = 'ANGLES'
      TKFRAME_-70109_RELATIVE   = 'DII_ITS'
      TKFRAME_-70109_ANGLES     = ( 0.0, 90.0, 60.0 )
      TKFRAME_-70109_AXES       = ( 3,    2,    1   )
      TKFRAME_-70109_UNITS      = 'DEGREES'

   \begintext


Comet Tempel body-fixed frame
--------------------------------------------------------

   The Tempel fixed frame is defined in the same as any other PCK frame:
   
      *  +Z along comet's North pole;
      
      *  +X along comet's prime meridian;
      
      *  +Y completes the right hand frame;
      
      *  the origin of this frame is at the center of the comet ellipsoid.
      
   As for any PCK frame orientation of this frame is computed by evaluating 
   corresponding rotation constants provided in a PCK file.

   \begindata

      FRAME_TEMPEL_FIXED     =  1000093
      FRAME_1000093_NAME     = 'TEMPEL_FIXED'
      FRAME_1000093_CLASS    =  2
      FRAME_1000093_CLASS_ID =  1000093
      FRAME_1000093_CENTER   =  1000093

      OBJECT_1000093_FRAME   = 'TEMPEL_FIXED'

   \begintext
