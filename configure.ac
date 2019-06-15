AC_INIT([gwgen], [0.1], [jed.kaplan@arve-research.com])

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/)

AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_FILES([src/Makefile])

AC_CONFIG_MACRO_DIR([m4])

AC_PREFIX_DEFAULT(["`pwd`"])

# Handle netcdf
AC_PATH_NETCDF

# For gwgen_grid, the netcdf V4 interface is REQUIRED.
if test $NETCDF_V4 = "no"; then
	echo "-------------------------------------------------------------------"
	echo "gwgen version 2.0 REQUIRES the netcdf library be version 4 or above,"
	echo "AND installed with HDF-5 support (i.e., the netcdf library must be"
	echo "compiled with the --enable-netcdf-4 flag). If you don't have netcdf-4"
	echo "with HDF-5 support, then please install it first before proceeding"
	echo "-------------------------------------------------------------------"
	exit -1
fi

#------------------------------------------------------------------------
# Check for the C and Fortran compilers. NOTE that this must come AFTER the netcdf check. 
# The netcdf libraries could have been compiled with, for example, Intel C and ifort;
# we have to specifically check if the same compiler used to compile the
# netcdf library exists.
#------------------------------------------------------------------------

AC_PROG_CC( $NETCDF_CC )
AC_PROG_FC( $NETCDF_FC )
AC_LANG([Fortran])
AC_FC_SRCEXT(f90)

AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
    [compile with MPI (parallelization) support. If none is found,
    MPI is not used. Default: auto])
],,[with_mpi=auto])

AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
  else
    AC_MSG_WARN([No MPI compiler found, won't use MPI.])
  fi
])

AX_F90_MODULE([netcdfmod], ["netcdf.mod"], [
  use netcdf
  status = nf90_close(ncid)
], [/usr/local/include:/sw/include], [$with_netcdfmod], [AC_MSG_ERROR(Can't find netcdf module.)])

save_FCFLAGS=$FCFLAGS
save_LDFLAGS=$LDFLAGS
FCFLAGS="$FCFLAGS $F90_MODULE_netcdfmod"

AX_F90_LIBRARY([netcdf], [libnetcdff*], [
  use netcdf
  status = nf90_close(ncid)
], [/usr/local/lib:/sw/lib], [$with_netcdf], [AC_MSG_ERROR(Can't find netcdf library.)])

LDFLAGS="$LDFLAGS $F90_LIBS_netcdf"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
