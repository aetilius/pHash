#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([pHash],[0.9.6],[support@phash.org])
AM_INIT_AUTOMAKE($PACKAGE_NAME,$PACKAGE_VERSION)
AC_CONFIG_SRCDIR([src/pHash.cpp])
AC_CONFIG_HEADERS([config.h pHash-config.h])

m4_include([m4/ax_pthread.m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

if test -z $DISTRO; then
   AC_CHECK_FILE(/etc/gentoo-release, [DISTRO="gentoo"])
   AC_CHECK_FILE(/etc/redhat-release, [DISTRO="redhat"])
   AC_CHECK_FILE(/etc/slackware-version, [DISTRO="slackware"])
   AC_CHECK_FILE(/etc/debian_version, [DISTRO="debian"])
   AC_CHECK_FILE(/etc/SuSErelease, [DISTRO="suse"])
fi


AC_ARG_ENABLE(pthread, [AS_HELP_STRING([--disable-pthread],
[pthread support @<:@default=no@:>@])])

AS_IF([test "x$enable_pthread" != "xno"],
[
	AX_PTHREAD([
	LIBS="$PTHREAD_LIBS $LIBS"
	CPPFLAGS="$CPPFLAGS $PTHREAD_CFLAGS"
	CC="$PTHREAD_CC"
	AC_DEFINE([HAVE_PTHREAD], [1], [configure with pthread support])
	])
])

AM_CONDITIONAL([HAVE_PTHREAD], [test "x$enable_pthread" != "xno"])

AC_ARG_ENABLE(debug, [AS_HELP_STRING([--enable-debug],
[compile with debugging support @<:@default=no@:>@])])

AS_IF([test "x$enable_debug" = "xyes"],
	[CXXFLAGS="`echo $CXXFLAGS' ' | sed 's/-O[0-9]//g'`"
	CXXFLAGS="$CXXFLAGS -O0 -g"
	LDFLAGS="`echo $LDFLAGS | sed 's/-g//g' | sed 's/-O[0-9]//g'`"],
	[CXXFLAGS="`echo $CXXFLAGS' -ffast-math' | sed 's/-g//g' | sed 's/-O[0-9]//g'`"
	CXXFLAGS="$CXXFLAGS -O3"
	CPPFLAGS="`echo $CPPFLAGS | sed 's/-g//g' | sed 's/-O[0-9]//g'`"
	LDFLAGS="`echo $LDFLAGS | sed 's/-g//g' | sed 's/-O[0-9]//g'`"])

AC_ARG_ENABLE(openmp, [AS_HELP_STRING([--enable-openmp],
[enable OpenMP support in pHash to use multiple cores/CPUs @<:@default=no@:>@])])

AS_IF([test "x$enable_openmp" = "xyes"],
[
	CPPFLAGS="$CPPFLAGS -Dcimg_use_openmp"
	if test $GCC = "yes"; then
		CXXFLAGS="$CXXFLAGS -fopenmp"
	fi
])

AC_ARG_ENABLE(java, [AS_HELP_STRING([--enable-java],
[compile java (JNI) binding for pHash @<:@default=no@:>@])])

AS_IF([test "x$enable_java" = "xyes"],
[
	AC_CHECK_HEADER([jni.h], [found_jni="y"], [found_jni="n"])

if test x"$found_jni" = x"n"; then
	AC_MSG_CHECKING([whether jni.h is in the current or src directory.])
	
	if [test `ls . src | grep "jni.h"` >/dev/null 2>&1]; then
		AC_MSG_RESULT([yes])
		found_java="y"
	else
		AC_MSG_RESULT([no])
		AC_MSG_CHECKING([if jni.h is in other system directories.])
		jni_header=`locate jni.h|grep include/jni.h|head -n1`
		
		if test x"$jni_header" != "x"; then
			AC_MSG_RESULT([yes])
			jni_header=`dirname $jni_header`
			jni_md_header=$jni_header"/"`uname|tr '[A-Z]' '[a-z]'`
			CPPFLAGS="$CPPFLAGS -I $jni_header -I $jni_md_header"
		else
			AC_MSG_RESULT([no])
		fi
	fi
	
fi

])

AM_CONDITIONAL([WITH_JAVA], [test "x$enable_java" = "xyes"])

AC_PROG_INSTALL
AC_LANG([C++])
AC_PROG_LN_S
AC_PROG_LIBTOOL

AC_SUBST([LIBTOOLS_DEPS])

LDFLAGS="$LDFLAGS -L/usr/local/lib"
CPPFLAGS="$CPPFLAGS -I/usr/local/include"

AC_DEFUN([AC_CHECK_CIMG], [
AC_CHECK_HEADER([CImg.h], [found_cimg="y"], [found_cimg="n"])

if test x"$found_cimg" = x"n"; then
	AC_MSG_CHECKING([whether CImg.h is in the current or src directory.])
	
	if [test `ls . src | grep "CImg.h"` >/dev/null 2>&1]; then
		AC_MSG_RESULT([yes])
		found_cimg="y"
	else
		AC_MSG_RESULT([no])
	fi
	
fi])
 
AC_DEFUN([AC_CHECK_FFMPEG],
[
AC_MSG_CHECKING([whether FFmpeg is present])
AC_CHECK_LIB([avcodec], [avcodec_alloc_frame], [], [AC_MSG_ERROR([

*** libavcodec not found.
You need FFmpeg. Get it at <http://ffmpeg.org/>])])

AC_CHECK_LIB([avutil], [av_log_set_level], [], [AC_MSG_ERROR([

*** libavutil not found.
You need FFmpeg. Get it at <http://ffmpeg.org/>])])

AC_CHECK_LIB([avformat], [av_read_frame], [], [AC_MSG_ERROR([

*** libavformat not found.
You need FFmpeg. Get it at <http://ffmpeg.org/>])])

AC_CHECK_LIB([swscale], [sws_getContext], [], [AC_MSG_ERROR([

*** libswscale not found.
You need FFmpeg. Get it at <http://ffmpeg.org/>])])
found_ffmpeg="y"])


AC_DEFUN([AC_CHECK_SNDFILE],
[
AC_MSG_CHECKING([whether sndfile lib is present])
AC_CHECK_LIB([sndfile], [sf_readf_float], [], [AC_MSG_ERROR([

*** libsndfile not found.
You need sndfile.  Get it at <http://www.mega-nerd.com/libsndfile/>])])
found_sndfile="y"])


AC_DEFUN([AC_CHECK_SAMPLERATE],
[
AC_MSG_CHECKING([whether samplerate lib is present])
AC_CHECK_LIB([samplerate],[src_process],[],[AC_MSG_ERROR([

*** libsamplerate not found.
You need samplerate. Get it at <http://www.mega-nerd.com/libsamplerate/>])])
found_samplerate="y"])


AC_DEFUN([AC_CHECK_MPG123],
[
AC_MSG_CHECKING([whether libmpg123 is present])
AC_CHECK_LIB([mpg123],[mpg123_init],[],[AC_MSG_ERROR([

*** libmpg123 not found.
You need libmpg123.
Get it at <http://www.mpg123.de/>])])
found_mpg123="y"])


AC_ARG_ENABLE([image-hash], [AS_HELP_STRING([--enable-image-hash],
[include support for perceptual image hashes @<:@default=yes@:>@])])

AS_IF([test x"$enable_image_hash" != "xno"],
[
	echo
	echo "*** Configuring image hash ***"
	echo
	AC_CHECK_CIMG()
	AC_DEFINE([HAVE_IMAGE_HASH], [1], [configure with image hash])
])

AM_CONDITIONAL(HAVE_IMAGE_HASH, [test "x$enable_image_hash" != "xno"])

AC_ARG_ENABLE([video-hash], [AS_HELP_STRING([--enable-video-hash],
[include support for perceptual video hashes @<:@default=yes@:>@])])

AS_IF([test "x$enable_video_hash" != "xno"],
[
	echo
	echo "*** Configuring video Hash ***"
	echo
	if test "$found_ffmpeg" != "y"; then
		AC_CHECK_FFMPEG() 

		case "$DISTRO" in
			debian) CPPFLAGS="$CPPFLAGS -I/usr/include/ffmpeg" ;;
			redhat) CPPFLAGS="$CPPFLAGS -I/usr/include/libavcodec -I/usr/include/libavformat -I/usr/include/libswscale" ;;
		esac
	fi
	if test "$found_cimg" != "y"; then
		AC_CHECK_CIMG()
	fi
	AC_DEFINE([HAVE_VIDEO_HASH],[1], [configure with video hash])
])

AM_CONDITIONAL(HAVE_VIDEO_HASH, [test "x$enable_video_hash" != "xno"])

AC_ARG_ENABLE([audio-hash], [AS_HELP_STRING([--enable-audio-hash],
[include support for perceptual audio hashes @<:@default=yes@:>@])])

AS_IF([test "x$enable_audio_hash" != "xno"],
	[
	echo
	echo "*** Configuring audio hash ***"
	echo
	if test "$found_sndfile" != "y"; then
		AC_CHECK_SNDFILE() 
	fi

	if test "$found_samplerate" != "y"; then
		AC_CHECK_SAMPLERATE()
	fi
	if test "$found_mpg123" != "y"; then
		AC_CHECK_MPG123()
	fi
	AC_DEFINE([HAVE_AUDIO_HASH],[1], [configure with audio hash])
	])

AM_CONDITIONAL(HAVE_AUDIO_HASH, [test "x$enable_audio_hash" != "xno"])

# Checks for libraries.
AC_ARG_WITH([libpng], AS_HELP_STRING([--without-libpng], [do not include libpng support]))
AC_ARG_WITH([libjpeg], AS_HELP_STRING([--without-libjpeg], [do not include libjpeg support]))
AC_ARG_WITH([fftw], AS_HELP_STRING([--without-fftw], [do not include FFTW support]))

AS_IF([test "x$enable_image_hash" != "xno" -o "x$enable_video_hash" != "xno"],
  AS_IF([test "x$with_libpng" != "xno"],
    AC_CHECK_LIB([png], [png_create_read_struct],
      [ CPPFLAGS="$CPPFLAGS -Dcimg_use_png" LIBS="$LIBS -lpng" ],
      AS_IF([test "x$with_libpng" = "xyes"], AC_MSG_ERROR([libpng support requested but library not found]))
    )
  )

  AS_IF([test "x$with_libjpeg" != "xno"],
    AC_CHECK_LIB([jpeg], [jpeg_read_header],
      [ CPPFLAGS="$CPPFLAGS -Dcimg_use_jpeg" LIBS="$LIBS -ljpeg" ],
      AS_IF([test "x$with_libjpeg" = "xyes"], AC_MSG_ERROR([libjpeg support requested but library not found]))
    )
  )
)

AS_IF([test "x$with_fftw" != "xno"],
  AC_CHECK_LIB([fftw3], [fftw_plan_dft_r2c_1d],
    [ CPPFLAGS="$CPPFLAGS -Dcimg_use_fftw3" LIBS="$LIBS -lfftw3" ],
    AS_IF([test "x$with_fftw" = "xyes"], AC_MSG_ERROR([FFTW support requested but library not found]))
  )
)

CFLAGS="$CXXFLAGS"

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_CHECK_TYPE([uint64_t],[],[AC_MSG_ERROR([

*** Required type uint64_t not found. 
pHash requires 64-bit integers for video and image hashes.])])

AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile pHash.pc bindings/Makefile bindings/java/Makefile])
AC_OUTPUT

