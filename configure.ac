#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([aligntk], [1.0], [aligntk@psc.edu], [http://mmbios.org/index.php/software/aligntk])
AC_CONFIG_SRCDIR([dt.h])
AC_CONFIG_HEADERS([config.h])

# Check for options
AC_ARG_WITH(x,
            [  --without-x             builds AlignTK without GUI programs],
            [], [with_x=yes])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_CHECK_PROG(MKDIR, [mkdir], [mkdir])
AC_CHECK_PROG(MPICC, [mpicc], [mpicc])

# Checks for libraries.

AC_LANG([C])
AC_CHECK_LIB([fftw3f], [fftwf_execute], TESTLIBS="$TESTLIBS -lfftw3f", AC_MSG_ERROR("Can't build without libfftw3f (does your system require a fftw3-devel package?)"), -lm)
AC_CHECK_LIB([jpeg], [jpeg_start_decompress], , AC_MSG_ERROR("Can't build without libjpeg (does your system require a libjpg-devel package?)"), -lm)
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([tiff], [TIFFOpen], , AC_MSG_ERROR("Can't build without libtiff (does your system require a libtiff-devel package?)"), -lm)
AC_CHECK_LIB([z], [gzopen], , AC_MSG_ERROR("Can't build without libz (does you system require a zlib-devel package?)"))

AC_LANG([C++])
AC_CHECK_LIB([fltk], [Fl::run])
AC_CHECK_LIB([fltk_gl], [Fl::gl_visual])

if test $with_x = "yes"; then
   AC_LANG([C])
   # X86_64 Linux machines have both 64 bit and 32 bit libraries. We need to choose the right set
   AC_CHECK_FILE([/usr/X11R6/lib64], [TESTLIBS="$TESTLIBS -L/usr/X11R6/lib64"], AC_CHECK_FILE([/usr/X11R6/lib], [TESTLIBS="$TESTLIBS -L/usr/X11R6/lib"]))
   AC_CHECK_LIB([X11], [XOpenDisplay], TESTLIBS="$TESTLIBS -lX11",, $TESTLIBS)
   AC_CHECK_LIB([GL], [glBegin], TESTLIBS="$TESTLIBS -lGL",, $TESTLIBS)
   AC_CHECK_LIB([GLU], [gluBeginPolygon], TESTLIBS="$TESTLIBS -lGLU",, $TESTLIBS)
   AC_SUBST(TARGETS, "nox_executables x_executables")
   AC_SUBST(INSTALL_TARGETS, "nox_install x_install")
else
   AC_SUBST(TARGETS, "nox_executables")
   AC_SUBST(INSTALL_TARGETS, "nox_install")
fi

# Checks for header files.
AC_CHECK_HEADERS([limits.h netdb.h netinet/in.h stdlib.h string.h sys/socket.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

AC_LANG([C++])
AC_CHECK_HEADERS([FL/Fl.H])
AC_CHECK_HEADERS([FL/fl_draw.H])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor gethostbyname gethostname gettimeofday localtime_r memmove memset mkdir sqrt strcasecmp strchr strdup strrchr])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
