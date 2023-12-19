#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([guinea-pig],[1.2.2],[barbara.dalena@cern.ch])

EXTERNAL_CFLAGS="$CFLAGS"
EXTERNAL_CXXFLAGS="$CXXFLAGS"

DEFAULT_CFLAGS="-O3"
DEFAULT_CXXFLAGS=""

AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([src/guineapigCPP.h])
AC_CONFIG_HEADERS([src/config.h])

if test "$prefix" = "NONE" ; then
  prefix="/usr/local"
fi

make_all=guinea

USE_FFTW2=false
USE_FFTW3=false
USE_ROOT=false

AC_ARG_ENABLE([fftw2], [AS_HELP_STRING([--enable-fftw2],[enable the fftw2 libraries for field computation (default is no)])], [USE_FFTW2=true], [USE_FFTW2=false])
AC_ARG_ENABLE([fftw3], [AS_HELP_STRING([--enable-fftw3],[enable the fftw3 libraries for field computation (default is no)])], [USE_FFTW3=true], [USE_FFTW3=false])
AC_ARG_ENABLE([root], [AS_HELP_STRING([--enable-root],[enable the ROOT libraries for outputs (default is no)])], [USE_ROOT=true], [USE_ROOT=false])


AC_ARG_WITH([fftwdir], [AS_HELP_STRING([--with-fftwdir=DIR],[FFTW home directory (such that `$DIR/lib' and  `$DIR/include' exist)])], [FFTW_HOME="$with_fftwdir"])
AC_ARG_WITH([rootdir], [AS_HELP_STRING([--with-rootdir=DIR],[ROOT home directory (such that `$DIR/lib' and  `$DIR/include' exist)])], [ROOTSYS="$with_rootdir"])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

CFLAGS="$DEFAULT_CFLAGS"
CXXFLAGS="$DEFAULT_CXXFLAGS"

#
# setting up hardware dependent options
#

CPPFLAGS="$CPPFLAGS -DUSE_NON_CONST"
CFLAGS="$CFLAGS -ffast-math -fno-keep-inline-functions -Wall -Wundef"

case "$target" in
i[[3456]]86-*-linux-*)
  LDFLAGS="$LDFLAGS"
  CFLAGS="$CFLAGS -fPIC"
  CPPFLAGS="$CPPFLAGS -malign-double"
;;
x86_64-*-linux-*)
  LDFLAGS="$LDFLAGS"
  CFLAGS="$CFLAGS -fPIC"
;;
*-apple-darwin*)
  CFLAGS="$CFLAGS"
  LDFLAGS="$LDFLAGS -L/sw/lib"
;;
esac

# Check for 64 bit
AC_CHECK_SIZEOF([void*])

if SIZEOF_VOID_P=8; then
   CFLAGS="$CFLAGS -m64 -DCOMPUTER64b"
fi

# Checks for libraries.
#
# checks for FFTW
#
if test "$USE_FFTW2" = "true" ; then

   test -z $FFTW_HOME  && AC_CHECK_PROG(FFTW_HOME, $FFTW_HOME/include/fftw.h)
     
     if test -z $FFTW_HOME; then
	AC_MSG_ERROR([\`$FFTW_HOME/include/fftw.h' not found. FFTW2 is required if to compile and run guinea. If you know that FFTW is properly installed, you can set the environment variable \`FFTW_HOME' to point to your installation directory.])
     fi

     FFTW_CFLAGS="-DUSE_FFTW2 -I$FFTW_HOME/include"
     FFTW_LDFLAGS="-L$FFTW_HOME/lib -lfftw -lm"
fi   

if test "$USE_FFTW3" = "true" ; then

   test -z $FFTW_HOME  && AC_CHECK_PROG(FFTW_HOME, $FFTW_HOME/include/fftw3.h)
     
     if test -z $FFTW_HOME; then
	AC_MSG_ERROR([\`$FFTW_HOME/include/fftw3.h' not found. FFTW3 is required if to compile and run guinea. If you know that FFTW is properly installed, you can set the environment variable \`FFTW_HOME' to point to your installation directory.])
     fi

     FFTW_CFLAGS="-DUSE_FFTW3 -I$FFTW_HOME/include"
     FFTW_LDFLAGS="-L$FFTW_HOME/lib -lfftw3 -lm"
fi   
#
# checks for ROOT
#

if test "$USE_ROOT" = "true" ; then

   test -z $ROOTSYS && AC_CHECK_PROG(ROOTSYS)
   
   if test -z $ROOTSYS; then
      AC_MSG_ERROR([\`ROOT' not found. ROOT is required to compile and run guinearoot. If you know that ROOT is properly installed, you can set the environment variable \`ROOTSYS' to point to your installation directory.])
   fi

   make_all ="$make_all guinearoot"
   ROOT_CFLAGS="-I$(ROOTSYS)/include"
   ROOT_LDFLAGS="-L$(ROOTSYS)/lib"

fi

# Checks for header files.
m4_warn([obsolete],
[The preprocessor macro `STDC_HEADERS' is obsolete.
  Except in unusual embedded environments, you can safely include all
  ISO C90 headers unconditionally.])dnl
# Autoupdate added the next two lines to ensure that your configure
# script's behavior did not change.  They are probably safe to remove.
AC_CHECK_INCLUDES_DEFAULT
AC_PROG_EGREP

AC_CHECK_HEADERS([stdlib.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow sqrt stpcpy strstr strtol])

#
# last operations before substitutions
#

if test -n "$EXTERNAL_CFLAGS" ; then
  CFLAGS="$CFLAGS $EXTERNAL_CFLAGS"
fi

if test -n "$EXTERNAL_CXXFLAGS" ; then
  CXXFLAGS="$CXXFLAGS $EXTERNAL_CXXFLAGS"
fi

if test "$USE_FFTW2" = "false" && test "$USE_FFTW3" = "false"; then
   FFTW_CFLAGS="-DUSE_FFT_LOCAL"   
fi

if test "$USE_ROOT" = "false" ; then
   ROOT_CFLAGS=""   
fi


AC_SUBST(make_all)
AC_SUBST(FFTW_CFLAGS)
AC_SUBST(ROOT_CFLAGS)
AC_SUBST(FFTW_LDFLAGS)
AC_SUBST(ROOT_LDFLAGS)

AC_CONFIG_FILES([Makefile src/Makefile])


AC_OUTPUT
