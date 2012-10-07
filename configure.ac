AC_INIT([allframe], [4.1.2], [pbs@nrc.dao.ca])

AC_CONFIG_SRCDIR([src/daophot.f])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])

AM_INIT_AUTOMAKE([foreign])

AC_PROG_F77
AC_ARG_ENABLE(debug,
   AC_HELP_STRING([--enable-debug], [enable memory debug mode (default=no)]),
   [enable_debug=$enableval],
   [enable_debug=no])
if test x"$enable_debug" = x"yes"; then
   AC_MSG_RESULT([yes])
   if test x$F77 = xg77; then
      FFLAGS="${FFLAGS} -fbounds-check"
      LDFLAGS="${LDFLAGS} --defsym,mem_=0"
   fi
fi

# automatically add -lm to all linking steps
AC_CHECK_LIB([m], [cos])

# Checks for cfitsio libraries
# TODO: check fortran bindings
AC_ARG_WITH(fits,
   AC_HELP_STRING([--without-fits], [disable FITS i/o capabilities]),
   [with_fits=$withval],
   [with_fits=yes])
if test x"$with_fits" = x"yes"; then
  AX_CHECK_PKG_LIB(
  	[cfitsio],
	[fitsio.h],
	[ffopen],,
	[AC_MSG_ERROR([Could not find cfitsio library])])
fi
AM_CONDITIONAL([WITH_FITS], [test x$with_fits = xyes])

# Checks for iraf libraries
AC_ARG_WITH(iraf,
   AC_HELP_STRING([--with-iraf], [enable IRAF i/o capabilities]),
   [with_iraf=$withval],
   [with_iraf=no])
if test x"$with_iraf" != x"no"; then
  AX_CHECK_PKG_LIB(
  	[iraf],
	[],
	[ftppre],,
	[AC_MSG_ERROR([Could not find iraf library])])
fi
AM_CONDITIONAL([WITH_IRAF], [test x$with_iraf = xyes])

AC_CONFIG_FILES([Makefile 
		src/Makefile])
AC_OUTPUT
