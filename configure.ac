AC_INIT([gwgen], [0.1], [jed.kaplan@arve-research.com])

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/)

AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_FILES([src/Makefile])

AC_CONFIG_MACRO_DIR([m4])

AC_PREFIX_DEFAULT(["`pwd`"])

AC_PROG_CC
AC_PROG_FC
AC_LANG([Fortran])
AC_FC_SRCEXT(f90)

AC_ARG_WITH(netcdf, [AS_HELP_STRING([--with-netcdf],
    [Specify the location of your netCDF installation. Default: /usr/local or /sw])
],,[with_netcdf=auto])

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
