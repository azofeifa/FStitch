#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.62])
AC_INIT([FStitch], [0.1], [joseph.azofeifa@colorado.edu])
AM_INIT_AUTOMAKE([foreign])
AM_SILENT_RULES([yes])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/viterbi.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11

# Checks for libraries.
AC_OPENMP
AC_SUBST(AM_CXXFLAGS,
  "$OPENMP_CXXFLAGS")

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
