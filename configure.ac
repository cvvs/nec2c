dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.62)
AC_INIT(nec2c, 0.9, mjtrangoni@gmail.com)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR([src/main.c])

AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([no-define])
AM_MAINTAINER_MODE
AC_PREFIX_DEFAULT([/usr/local])
test "x$prefix" = xNONE && prefix="/usr/local"

dnl Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_CHECK_TOOL(AR, ar, :)

dnl Check math library
AC_CHECK_LIB(m,main)

dnl Checks for library functions.
AC_CHECK_FUNCS([floor pow exp sqrt log10 log cos sin acos asin tan atan sinh \
		cosh tanh fabs modf atan2 jn yn erf erfc round trunc acosh \
		asinh strdup strerror strchr])

dnl Checks for complex classes and functions.
AC_CHECK_HEADERS([complex.h stdio.h signal.h math.h stdlib.h unistd.h \
		  string.h fcntl.h errno.h time.h sys/types.h sys/times.h \
		  omp.h],
	 	[mypj_found_int_headers=yes; break;])

AS_IF([test "x$mypj_found_int_headers" != "xyes"],
      	[AC_MSG_ERROR([Unable to find the standard headers to nec2c.h])])


dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

dnl Check for type sizes.

AC_CHECK_SIZEOF([short])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])

type32="int"
if test $ac_cv_sizeof_short -eq 4; then
	type32=short
elif test $ac_cv_sizeof_int -eq 4; then
	type32=int
elif test $ac_cv_sizeof_long -eq 4; then
      type32=long
	fi

dnl Check for double type.

AC_ARG_ENABLE([double],
AS_HELP_STRING([--enable-double=TYPE],
	[type of double representation @<:@default=double@:>@]),
		[case "$enableval" in
       		no | yes | double)
     		double_type="double"
     		double_size=8
     		;;
   		ldouble)
double_type="long double"
AC_CHECK_SIZEOF([long double],12)
double_size=$ac_cv_sizeof_long_double
;;
float)
double_type="float"
double_size=4
;;
*)
double_type="none"
double_size=0
;;
esac],
[double_type="double"; double_size=8])
if test "$double_type" = none; then
       AC_MSG_ERROR([valid types are double, float and ldouble.])
fi       

dnl Check OpenMP support
AC_OPENMP

AC_OUTPUT(Makefile src/Makefile)

dnl Print results.
AC_MSG_RESULT([])
AC_MSG_RESULT([ $PACKAGE version $VERSION configured successfully.])
AC_MSG_RESULT([])

