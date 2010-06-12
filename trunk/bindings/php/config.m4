dnl
dnl $ Id: $
dnl

PHP_ARG_WITH(pHash, whether pHash is available,[  --with-pHash[=DIR]   With pHash support])


if test "$PHP_PHASH" != "no"; then
  PHP_REQUIRE_CXX
  AC_LANG_CPLUSPLUS
  PHP_ADD_LIBRARY(stdc++,,PHASH_SHARED_LIBADD)


  if test -r "$PHP_PHASH/include/pHash.h"; then
	PHP_PHASH_DIR="$PHP_PHASH"
  else
	AC_MSG_CHECKING(for pHash in default path)
	for i in /usr /usr/local; do
	  if test -r "$i/include/pHash.h"; then
		PHP_PHASH_DIR=$i
		AC_MSG_RESULT(found in $i)
		break
	  fi
	done
	if test "x" = "x$PHP_PHASH_DIR"; then
	  AC_MSG_ERROR(not found)
	fi
  fi

  PHP_ADD_INCLUDE($PHP_PHASH_DIR/include)

  export OLD_CPPFLAGS="$CPPFLAGS"
  export CPPFLAGS="$CPPFLAGS $INCLUDES -DHAVE_PHASH"
  AC_CHECK_HEADER([pHash.h], [], AC_MSG_ERROR('pHash.h' header not found))
  AC_CHECK_HEADER([pHash-config.h], [], AC_MSG_ERROR('pHash-config.h' header not found))
  AC_CHECK_HEADER([audiophash.h], [], AC_MSG_ERROR('audiophash.h' header not found))
  AC_CHECK_HEADER([CImg.h], [], AC_MSG_ERROR('CImg.h' header not found))
  PHP_SUBST(PHASH_SHARED_LIBADD)


  PHP_CHECK_LIBRARY(avcodec, avcodec_alloc_frame,
  [
	PHP_ADD_LIBRARY_WITH_PATH(avcodec, $PHP_PHASH_DIR/lib, PHASH_SHARED_LIBADD)
  ],[
	AC_MSG_ERROR([wrong avcodec lib version or lib not found])
  ],[
	-L$PHP_PHASH_DIR/lib
  ])

  PHP_CHECK_LIBRARY(avutil, av_log_set_level,
  [
	PHP_ADD_LIBRARY_WITH_PATH(avutil, $PHP_PHASH_DIR/lib, PHASH_SHARED_LIBADD)
  ],[
	AC_MSG_ERROR([wrong avutil lib version or lib not found])
  ],[
	-L$PHP_PHASH_DIR/lib
  ])

  PHP_CHECK_LIBRARY(avformat, av_read_frame,
  [
	PHP_ADD_LIBRARY_WITH_PATH(avformat, $PHP_PHASH_DIR/lib, PHASH_SHARED_LIBADD)
  ],[
	AC_MSG_ERROR([wrong avformat lib version or lib not found])
  ],[
	-L$PHP_PHASH_DIR/lib
  ])

  PHP_CHECK_LIBRARY(swscale, sws_getContext,
  [
	PHP_ADD_LIBRARY_WITH_PATH(swscale, $PHP_PHASH_DIR/lib, PHASH_SHARED_LIBADD)
  ],[
	AC_MSG_ERROR([wrong swscale lib version or lib not found])
  ],[
	-L$PHP_PHASH_DIR/lib
  ])

  PHP_CHECK_LIBRARY(sndfile, sf_readf_float,
  [
	PHP_ADD_LIBRARY_WITH_PATH(sndfile, $PHP_PHASH_DIR/lib, PHASH_SHARED_LIBADD)
  ],[
	AC_MSG_ERROR([wrong sndfile lib version or lib not found])
  ],[
	-L$PHP_PHASH_DIR/lib
  ])

  PHP_CHECK_LIBRARY(samplerate, src_process,
  [
	PHP_ADD_LIBRARY_WITH_PATH(samplerate, $PHP_PHASH_DIR/lib, PHASH_SHARED_LIBADD)
  ],[
	AC_MSG_ERROR([wrong samplerate lib version or lib not found])
  ],[
	-L$PHP_PHASH_DIR/lib
  ])
  export CPPFLAGS="$OLD_CPPFLAGS"

  export OLD_CPPFLAGS="$CPPFLAGS"
  export CPPFLAGS="$CPPFLAGS $INCLUDES -DHAVE_PHASH"

  AC_MSG_CHECKING(PHP version)
  AC_TRY_COMPILE([#include <php_version.h>], [
#if PHP_VERSION_ID < 40000
#error  this extension requires at least PHP version 4.0.0
#endif
],
[AC_MSG_RESULT(ok)],
[AC_MSG_ERROR([need at least PHP 4.0.0])])

  export CPPFLAGS="$OLD_CPPFLAGS"


  PHP_SUBST(PHASH_SHARED_LIBADD)
  AC_DEFINE(HAVE_PHASH, 1, [ ])

  PHP_NEW_EXTENSION(pHash, pHash.cpp , $ext_shared)

fi

