AC_INIT(OTCETERA, 0.0.01-dev, mtholder@gmail.com)
AC_CONFIG_HEADERS([config.h])
AC_PREREQ(2.57)
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR(config)
AM_INIT_AUTOMAKE([1.6 foreign dist-zip tar-ustar filename-length-max=299])
AC_CONFIG_SRCDIR([ChangeLog])
# Checks for programs.
AC_PROG_LN_S
AC_PROG_INSTALL
AC_PROG_LIBTOOL

################################################################################
# calling AC_PROG_CXX resets CXXFLAGS, we use our own flags set in the
#	the AX_CXXFLAGS macro below.
# So we cache the argument to configure
#	here in ARG_CXX_FLAGS (so that our macro does not override them)
################################################################################
ARG_CXX_FLAGS="$CXXFLAGS"
ARG_CPP_FLAGS="$CPPFLAGS"
AC_PROG_CXX
CXXFLAGS="$ARG_CXX_FLAGS"
AX_CXXFLAGS
AC_MSG_CHECKING([whether C++ compiler works])
AC_LANG(C++)
AC_LANG_PUSH(C++)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <vector>
]], [[]])], [compiler_works='yes'], [compiler_works='no'])
AC_MSG_RESULT([$compiler_works])
AC_LANG_POP()
if test x$compiler_works = xno
then
    AC_MSG_ERROR([No functioning C++ compiler found])
fi
AC_SUBST(CXXFLAGS)
CPPFLAGS="-I\$(top_srcdir) $CPPFLAGS $ARG_CPP_FLAGS"
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])
# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_STRUCT_TM
AC_FUNC_MKTIME
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memmove strchr strtol])
BOOST_REQUIRE([1.54])
################################################################################
# the install-check requires python with subprocess to actually perform a check
# these tests should not cause failure of configure if python is not found.
################################################################################
AC_PATH_PROG([PYTHON], [python], [], [$PATH])
AC_PYTHON_MODULE([subprocess], [], [])
AC_SUBST(PYTHON)

if test -z "$PYTHON"
then
	MAYBE_TEST_DIR=
else
	MAYBE_TEST_DIR="test"
fi
AC_SUBST(MAYBE_TEST_DIR)
AC_PROG_MAKE_SET
AC_CONFIG_FILES([					\
	Makefile 						\
	otc/Makefile					\
	test/Makefile \
	tools/Makefile \
	otceterav0.0.pc \
	])
AC_OUTPUT
