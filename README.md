This program interpolates fields from an MPAS, FV3, or lat-lon mesh (the "source" mesh)
to a regular WRF grid (the "target" grid), as specified by a geogrid.exe file.

Right now, the "target" grid file MUST be a file as run through geogrid.exe.

Interpolation method is specified by file (see 'example' directory)

Compiling instructions: 

1) Edit build.sh script.  Needed libraries: NETCDF, ESMF, MPI.
2) Run build.sh

Code is based on Larissa Reames's "MPASSIT" program, but is expanded for other input
  data aside from MPAS, and the output is simplified.  Portions of the code are
  also different that allow for somewhat greater user run-time flexibility.
