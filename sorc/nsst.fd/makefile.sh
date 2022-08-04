#!/bin/sh
export FCOMP=ifort

export FFLAGS="-O3 -fp-model precise -g -r8 -i4"
export NETCDF_INCLUDE="-I${NETCDF}/include"
export NETCDF_LDFLAGS_F="-L${NETCDF}/lib -lnetcdf -lnetcdff"

make build
make install
make clean

err=$?
if [ $err -ne 0 ]; then
  echo ERROR BUILDING nst_to_rtg      
  exit 2
fi

exit
