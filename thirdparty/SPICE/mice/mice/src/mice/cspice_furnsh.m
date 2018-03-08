%-Abstract
%
%   CSPICE_FURNSH loads SPICE kernel files into MATLAB.
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
%      file   the name of a SPICE kernel file. The file may be either binary
%             or text. If the file is a binary SPICE kernel it will
%             be loaded into the appropriate SPICE subsystem.  If
%             'file' is a SPICE text kernel it will be loaded into the
%             kernel pool. If 'file' is a SPICE meta-kernel containing
%             initialization instructions (through use of the
%             correct kernel pool variables), the files specified in
%             those variables will be loaded into the appropriate
%             SPICE subsystem.
%
%             The SPICE text kernel format supports association of
%             names and data values using a "keyword = value"
%             format. The keyword-value pairs thus defined are
%             called "kernel variables."
%
%             While any information can be placed in a text kernel
%             file, the following string valued kernel variables are
%             recognized by SPICE as meta-kernel keywords:
%
%                  KERNELS_TO_LOAD
%                  PATH_SYMBOLS
%                  PATH_VALUES
%
%             Each kernel variable is discussed below.
%
%             KERNELS_TO_LOAD   is a list of SPICE kernels to be
%                               loaded into a program.  If file
%                               names do not fit within the kernel
%                               pool 80 character limit, they may be
%                               continued to subsequent array
%                               elements by placing the continuation
%                               character ('+') at the end of an
%                               element and then placing the
%                               remainder of the file name in the
%                               next array element.  (See the
%                               examples below for an illustration
%                               of this technique or consult the
%                               routine STPOOL for further details.)
%
%                               You may use one or more PATH_SYMBOL
%                               assignments (see below) to specify
%                               strings to be substituted for some
%                               part of a file name.
%
%             PATH_SYMBOLS      is a list of strings (without
%                               embedded blanks) which if
%                               encountered following the '$'
%                               character will be replaced with the
%                               corresponding PATH_VALUES string.
%                               Note that PATH_SYMBOLS are
%                               interpreted only in values
%                               associated with the KERNELS_TO_LOAD
%                               variable. There must be a one-to-one
%                               correspondence between the values
%                               supplied for PATH_SYMBOLS and
%                               PATH_VALUES. For the purpose of
%                               determining this correspondence, any
%                               path value that is continued over
%                               multiple array elements counts as a
%                               single value.
%
%             PATH_VALUES       is a list of expansions to use when
%                               PATH_SYMBOLS are encountered. If
%                               path values do not fit within the
%                               kernel pool 80 character limit, they
%                               may be continued in the same way as
%                               file names (see the KERNELS_TO_LOAD
%                               description above).
%
%             These kernel pool variables persist within the kernel
%             pool only until all kernels associated with the
%             variable KERNELS_TO_LOAD have been loaded.  Once all
%             specified kernels have been loaded, the variables
%             KERNELS_TO_LOAD, PATH_SYMBOLS and PATH_VALUES are
%             removed from the kernel pool.
%
%             [n,m] = size(file); char = class(file)
%
%                or
%
%             [1,n] = size(file); cell = class(file)
%
%   the call:
%
%      cspice_furnsh( file )
%
%   loads 'file' into the SPICE kernel system.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%
%     Example 1
%     ---------
%
%     Load the leapseconds kernel naif0010.tls and the planetary
%     ephemeris SPK file de405s.bsp.
%
%        cspice_furnsh( 'naif0010.tls' )
%        cspice_furnsh( 'de405s.bsp' )
%
%     Load a vector of kernel names or a cell of kernel names.
%
%        cspice_furnsh( { 'naif0008.tls', 'de405s.bsp' } )
%
%     Example 2
%     ---------
%
%     This example illustrates how you could create a meta-kernel file
%     for a program that requires several text and binary kernels.
%
%     First create a list of the kernels you need in a text file as
%     shown below.
%
%        \begintext
%
%           Here are the SPICE kernels required for my application
%           program.
%
%           Note that kernels are loaded in the order listed. Thus we
%           need to list the highest priority kernel last.
%
%
%        \begindata
%
%        KERNELS_TO_LOAD = (
%
%              '/home/mydir/kernels/spk/lowest_priority.bsp',
%              '/home/mydir/kernels/spk/next_priority.bsp',
%              '/home/mydir/kernels/spk/highest_priority.bsp',
%              '/home/mydir/kernels/text/leapsecond.ker',
%              '/home/mydir/kernels+',
%              '/custom+',
%              '/kernel_data/constants.ker',
%              '/home/mydir/kernels/text/sclk.tsc',
%              '/home/mydir/kernels/ck/c-kernel.bc' )
%
%
%     Note that the file name
%
%        /home/mydir/kernels/custom/kernel_data/constants.ker
%
%     is continued across several lines in the right hand side of the
%     assignment of the kernel variable KERNELS_TO_LOAD.
%
%     Once you've created your list of kernels, call FURNSH near the
%     beginning of your application program to load the meta-kernel
%     automatically at program start up.
%
%        cspice_furnsh( 'meta.tm' )
%
%     This will cause each of the kernels listed in your meta-kernel
%     to be loaded.
%
%
%     Example 3
%     ---------
%
%     This example illustrates how you can simplify the previous
%     kernel list by using PATH_SYMBOLS.
%
%
%        \begintext
%
%           Here are the SPICE kernels required for my application
%           program.
%
%           We are going to let A substitute for the directory that
%           contains SPK files; B substitute for the directory that
%           contains C-kernels; and C substitute for the directory that
%           contains text kernels.  And we'll let D substitute for
%           a "custom" directory that contains a special planetary
%           constants kernel made just for our mission.
%
%           Note that our PATH_VALUES and the corresponding
%           PATH_SYMBOLS must be listed in the same order.
%
%
%        \begindata
%
%        PATH_VALUES  = ( '/home/mydir/kernels/spk',
%                         '/home/mydir/kernels/ck',
%                         '/home/mydir/kernels/text',
%                         '/home/mydir/kernels/custom/kernel_data' )
%
%        PATH_SYMBOLS = ( 'A',
%                         'B',
%                         'C',
%                         'D'  )
%
%        KERNELS_TO_LOAD = (  '$A/lowest_priority.bsp',
%                             '$A/next_priority.bsp',
%                             '$A/highest_priority.bsp',
%                             '$C/leapsecond.ker',
%                             '$D/constants.ker',
%                             '$C/sclk.tsc',
%                             '$B/c-kernel.bc'         )
%
%
%     Example 4
%     ---------
%
%     This example illustrates continuation of path values. The
%     meta-kernel shown here is a modified version of that from
%     example 3.
%
%        \begintext
%
%           Here are the SPICE kernels required for my application
%           program.
%
%           We are going to let A substitute for the directory that
%           contains SPK files; B substitute for the directory that
%           contains C-kernels; and C substitute for the directory that
%           contains text kernels.  And we'll let D substitute for
%           a "custom" directory that contains a special planetary
%           constants kernel made just for our mission.
%
%           Note that our PATH_VALUES and the corresponding
%           PATH_SYMBOLS must be listed in the same order.
%
%           The values for path symbols A and D are continued over
%           multiple lines.
%
%        \begindata
%
%        PATH_VALUES  = ( '/very_long_top_level_path_name/mydir/+',
%                         'kernels/spk',
%                         '/home/mydir/kernels/ck',
%                         '/home/mydir/kernels/text',
%                         '/very_long_top_level_path_name+',
%                         '/mydir/kernels/custom+',
%                         '/kernel_data'                )
%
%        PATH_SYMBOLS = ( 'A',
%                         'B',
%                         'C',
%                         'D'  )
%
%        KERNELS_TO_LOAD = (  '$A/lowest_priority.bsp',
%                             '$A/next_priority.bsp',
%                             '$A/highest_priority.bsp',
%                             '$C/leapsecond.ker',
%                             '$D/constants.ker',
%                             '$C/sclk.tsc',
%                             '$B/c-kernel.bc'         )
%
%-Particulars
%
%   Text kernels input to this routine need not have native line
%   terminators for the platform. Lower level CSPICE routines can
%   read and process non-native text files. This functionality does
%   not exist in the FORTRAN SPICELIB.
%
%   Kernel pool variable names are restricted to a length of 32
%   characters or less.
%
%   In this version of the toolkit the maximum number of kernels that
%   can loaded together is limited to 5300. Each time a kernel is loaded
%   via cspice_furnsh, an internal kernel database entry is created for
%   that kernel. If a meta-kernel is loaded, a database entry is created
%   for the meta-kernel itself and for all files referenced in the
%   meta-kernel's KERNELS_TO_LOAD specification. Unloading a kernel or
%   meta-kernel deletes database entries created when the file was
%   loaded.
% 
%   The value above is an upper bound on number of SPICE kernels that
%   can be loaded at any time via cspice_furnsh, but the number of
%   kernels that can be loaded may be smaller, since re-loading a loaded
%   kernel or meta-kernel results in creation of additional database
%   entries.
% 
%   Kernels loaded via cspice_furnsh are subject to constraints imposed
%   by lower-level subsystems. The binary kernel systems (SPK, CK,
%   binary PCK, EK, and DSK) have their own limits on the maximum number
%   of kernels that may be loaded.
% 
%   The total number of DAF-based files (this set includes SPKs, CKs,
%   and binary PCKs) and DAS-based files (this set includes EKs and
%   DSKs) that may be loaded at any time may not exceed 5000. This limit
%   applies whether the files are loaded via cspice_furnsh or
%   lower-level loaders such as cspice_dafopr. File access performance
%   normally will degrade slightly as the number of loaded kernels
%   increases.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine furnsh_c.
%
%   MICE.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.0.3, 01-FEB-2017, BVS (JPL)
%
%      Added mention of the maximum number of kernels that can be loaded 
%      together. 
%
%   -Mice Version 1.0.2, 03-DEC-2013, EDW (JPL), SCK (JPL)
%
%      Expanded description of 'file' argument to match
%      all other SPICE language descriptions.
%
%      I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.1, 10-FEB-2010, EDW (JPL)
%
%      Added mention of the length restriction on the kernel pool variable
%      names.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   Load SPICE data from a list of items
%
%-&

function cspice_furnsh(file)

   switch nargin
      case 1

         file = zzmice_str(file);

      otherwise

         error ( 'Usage: cspice_furnsh(_`file`_)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('furnsh_c',file);
   catch
      rethrow(lasterror)
   end


