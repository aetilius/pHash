#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([pHash], [0.0.5], [support@phash.org])
AM_INIT_AUTOMAKE($PACKAGE_NAME,$PACKAGE_VERSION)
AC_CONFIG_SRCDIR([test_audiophash_main.cpp])
AC_CONFIG_HEADERS([config.h])
# Checks for programs.

AC_PROG_CXX
AC_PROG_CC

AC_MSG_CHECKING(whether to enable debugging)
AC_ARG_ENABLE(debug, [  --enable-debug[=no/yes] Turn on debugging [default=no]],
[ if test "x$enableval" = "xyes"; then
	CXXFLAGS="-g"
  else
	CXXFLAGS="-O2 -ffast-math"
  fi
])
AC_PROG_INSTALL
AC_PROG_LN_S
m4_defun([_LT_AC_LANG_F77_CONFIG], [:])
AC_PROG_LIBTOOL


# Checks for libraries.
AC_CHECK_LIB([avcodec], [avcodec_alloc_frame])
AC_CHECK_LIB([avformat], [avformat_alloc_context])
AC_CHECK_LIB([avutil], [avutil_version])
AC_CHECK_LIB([fftw3], [fftw_plan_dft_r2c_1d])
AC_CHECK_LIB([jpeg], [jpeg_read_header])
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([png], [png_create_read_struct])
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([swscale], [sws_getContext])

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor gettimeofday memmove memset pow sqrt strcasecmp strdup strncasecmp])
AC_CONFIG_FILES([Makefile pHash.pc])
AC_OUTPUT
