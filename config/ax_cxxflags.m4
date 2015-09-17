dnl @synopsis AX_CXXFLAGS()
dnl @author Mark Holder based on examples by Dave Swofford
dnl @version 2006-01-02
dnl @license AllPermissive
AC_DEFUN([AX_CXXFLAGS],
[
	dnl when using libtool to build shared libraries, we should use the
	dnl 	fPIC flag with g++.  MTH added, and I am not sure how many 
	dnl		compilers use this flag.  We will probably need to make this
	dnl 	conditional on libtool and compiler
	dnl CXXFLAGS_SHARED_LIB="-fPIC"
	dnl		COMMENTED out.  Looking at libtool script, I realized that this
	dnl		flag is supposed to be a libtool internal flag.  So I should figure
	dnl 	out why the libtool configure is not working correctly instead of
	dnl		adding the flag to all C++ compilation.
	
	dnl Check for debugging mode.
	AC_ARG_ENABLE(
		debugging,
		AC_HELP_STRING(
			[--enable-debugging],
			[build for debugging]
			), 
		, 
		[enable_debugging=no]
		)
	if test "$enable_debugging" = yes
	then
		AC_MSG_NOTICE([*** NOTE: debugging is enabled; optimization is suppressed!])
	fi

	
	if test "$enable_debugging" = yes; then
		CXXFLAGS_OPTIM_SPEED="-O0"
		CXXFLAGS="$CXXFLAGS -g"
	else
		CXXFLAGS_OPTIM_SPEED="-O"
	fi

	if test "$CC" = "icc"; then
			#	Intel C compiler for Linux
		if test "$enable_debugging" = no; then
			CXXFLAGS_OPTIM_SPEED="-O3"
			CXXFLAGS_OPTIM_SIZE="-O2"
		fi
	elif test "$CC" = "ccc"; then
			#	Compaq C compiler for Linux
		if test "x$arch" = "x"; then
			arch="host"
		fi
		if test "$enable_debugging" = no; then
			CXXFLAGS_OPTIM_SPEED="-fast -inline speed -arch $arch"
			CXXFLAGS_OPTIM_SIZE="-fast -inline size -unroll 1 -arch $arch"
		fi
	elif test "$CC" = "xlc"; then
			#	IBM XL C compiler
		if test "x$arch" = "x"; then
			arch="auto"
		fi
		if test "$enable_debugging" = no; then
			CXXFLAGS_OPTIM_SPEED="-O3 -qarch=$arch"
			CXXFLAGS_OPTIM_SIZE="-O3 -qarch=$arch"
		fi
	elif test "$GCC" = "yes" ; then
		if test "$enable_debugging" = yes; then
				#	Suppress warnings about possibly uninitialized variables but show everything else (used for
				#   development, but these warnings should also not trip for release builds)
			CXXFLAGS_WARNINGS="$CXXFLAGS_WARNINGS -Wall -Wno-uninitialized -Wimplicit -Wreturn-type -Wunused -Wredundant-decls -Wcast-align -Wcomment -Wextra"
		else
				#	Just suppress warnings about possibly uninitialized variables	
			CXXFLAGS_WARNINGS="$CXXFLAGS_WARNINGS -Wno-uninitialized"
			CXXFLAGS_OPTIM_SPEED="-O3 -ffast-math"
			CXXFLAGS_OPTIM_SIZE="-Os -ffast-math"
		fi
		case "$build_os" in
			darwin*) CXXFLAGS_WARNINGS="$CXXFLAGS_WARNINGS -Wno-long-double";;
			*);;
		esac	
	fi
	if test "x$CXXFLAGS_OPTIM_SIZE" = "x"; then
		CXXFLAGS_OPTIM_SIZE=$CXXFLAGS_OPTIM_SPEED
	fi
	
	CXXFLAGS="$CXXFLAGS $CXXFLAGS_OPTIM_SPEED $CXXFLAGS_WARNINGS $CXXFLAGS_SHARED_LIB"

])
