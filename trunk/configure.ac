#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.63)
AC_INIT([pHash],[0.5.0],[support@phash.org])
AM_INIT_AUTOMAKE($PACKAGE_NAME,$PACKAGE_VERSION)
AC_CONFIG_SRCDIR([src/pHash.cpp])
AC_CONFIG_HEADERS([config.h])
# Checks for programs.

AC_PROG_CXX
AC_PROG_CC
AM_PROG_GCJ

AC_ARG_ENABLE(debug, [AS_HELP_STRING([--enable-debug],
[compile with debugging support @<:@default=no@:>@])], 
debug=$enableval, debug=no)
AS_IF([test x"$debug" != x"no"],
	[CXXFLAGS="`echo $CXXFLAGS' ' | sed 's/-O[0-9]//g'`"
	CXXFLAGS="$CXXFLAGS -O0 -g"],
	[CXXFLAGS="`echo $CXXFLAGS' -ffast-math' | sed 's/-g//g' | sed 's/-O[0-9]//g'`"
	CXXFLAGS="$CXXFLAGS -O3"])

AC_ARG_ENABLE(java, [AS_HELP_STRING([--enable-java],
[compile java (JNI) binding for pHash @<:@default=no@:>@])], 
java=yes, java=no)
AM_CONDITIONAL([WITH_JAVA], test x$java != xno)
AC_PROG_INSTALL
AC_PROG_LN_S
m4_defun([_LT_AC_LANG_F77_CONFIG], [:])

AC_PROG_LIBTOOL
LT_INIT
AC_SUBST([LIBTOOLS_DEPS])

LDFLAGS="$LDFLAGS -L/usr/local/lib"
CPPFLAGS="$CPPFLAGS -I/usr/local/include"
AC_DEFUN([AC_CHECK_CIMG], [
AC_MSG_CHECKING([whether CImg.h is present in system directories.])
AC_MSG_RESULT([])
AC_CHECK_HEADER([CImg.h], [found_cimg="y"], [found_cimg="n"])

if test "$found_cimg" = "n"; then
	AC_MSG_CHECKING([whether CImg.h is in the src directory.])
	
	if [test `ls src | grep "CImg.h"` >/dev/null 2>&1]; then
		AC_MSG_RESULT([yes])
		found_cimg="y"
	else
		AC_MSG_RESULT([no])
	fi
	
fi])
AC_DEFUN([AC_CHECK_FFTW3],
[

AC_MSG_CHECKING([whether FFTW3 is present])
AC_CHECK_LIB([fftw3], [fftw_plan_dft_r2c_1d], [], [AC_MSG_ERROR([

*** libfftw3 not found.
You need FFTW3. Get it at <http://fftw.org/>])])
found_fftw="y"])

AC_DEFUN([AC_CHECK_FFMPEG],
[
AC_MSG_CHECKING([whether FFmpeg is present])
AC_CHECK_LIB([avcodec], [avcodec_alloc_frame], [], [AC_MSG_ERROR([

*** libavcodec not found.
You need FFmpeg. Get it at <http://ffmpeg.org/>])])

AC_CHECK_LIB([avformat], [av_read_frame], [], [AC_MSG_ERROR([

*** libavformat not found.
You need FFmpeg. Get it at <http://ffmpeg.org/>])])

AC_CHECK_LIB([swscale], [sws_getContext], [], [AC_MSG_ERROR([

*** libswscale not found.
You need FFmpeg. Get it at <http://ffmpeg.org/>])])
found_ffmpeg="y"])

AC_ARG_ENABLE([image-hash], [AS_HELP_STRING([--enable-image-hash],[include support for perceptual image hashes @<:@default=yes@:>@])],
[image_hash=$enableval], [image_hash=yes])
AS_IF([test x"$image_hash" != x"no"],
	[
	echo
	echo "*** Configuring image hash ***"
	echo
	AC_CHECK_CIMG()
	CPPFLAGS="$CPPFLAGS -DHAVE_IMAGE_HASH=1"])
AM_CONDITIONAL(HAVE_IMAGE_HASH, test x$image_hash != xno)

AC_ARG_ENABLE([video-hash], [AS_HELP_STRING([--enable-video-hash],[include support for perceptual video hashes @<:@default=yes@:>@])],
[video_hash=$enableval], [video_hash=yes])
AS_IF([test x"$video_hash" != x"no"],
	[
	echo
	echo "*** Configuring video Hash ***"
	echo
	if test "$found_ffmpeg" != "y"; then
		AC_CHECK_FFMPEG() 
	fi
	if test "$found_cimg" != "y"; then
		AC_CHECK_CIMG()
	fi
	CPPFLAGS="$CPPFLAGS -DHAVE_VIDEO_HASH=1"])
AM_CONDITIONAL(HAVE_VIDEO_HASH, test x$video_hash != xno)

AC_ARG_ENABLE([audio-hash], [AS_HELP_STRING([--enable-audio-hash],[include support for perceptual audio hashes @<:@default=yes@:>@])],
[audio_hash=$enableval], [audio_hash=yes])
AS_IF([test x"$audio_hash" != x"no"],
	[
	echo
	echo "*** Configuring audio hash ***"
	echo
	if test "$found_ffmpeg" != "y"; then
		AC_CHECK_FFMPEG() 
	fi
	if test "$found_fftw" != "y"; then
		AC_CHECK_FFTW3() 
	fi
	CPPFLAGS="$CPPFLAGS -DHAVE_AUDIO_HASH=1"])
AM_CONDITIONAL(HAVE_AUDIO_HASH, test x$audio_hash != xno)

# Checks for libraries.


AC_CHECK_LIB([jpeg], [jpeg_read_header], [use_jpeg="y"])
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([png], [png_create_read_struct], [use_png="y"])
AC_CHECK_LIB([pthread], [pthread_create])

AS_IF([test "$use_jpeg" = "y"],
[
	CPPFLAGS="$CPPFLAGS -Dcimg_use_jpeg"
	LDFLAGS="$LDFLAGS -ljpeg"])
AS_IF([test "$use_png" = "y"],
[
	CPPFLAGS="$CPPFLAGS -Dcimg_use_png"
	LDFLAGS="$LDFLAGS -lpng"])


# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_CHECK_TYPE([uint64_t],[],[AC_MSG_ERROR([

*** Required type uint64_t not found. 
pHash requires 64-bit integers for video and image hashes.])])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor gettimeofday memmove memset pow sqrt strcasecmp strdup strncasecmp])
AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile pHash.pc bindings/Makefile bindings/java/Makefile])
AC_OUTPUT
