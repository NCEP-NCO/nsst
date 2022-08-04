#! /usr/bin/env bash
set -eux

module reset
source ../versions/build.ver

module load PrgEnv-intel/$PrgEnv_intel_ver
module load intel/$intel_ver
module load craype/$craype_ver
module load cray-mpich/$cray_mpich_ver
module load w3nco/$w3nco_ver
module load bacio/$bacio_ver
module load nemsio/$nemsio_ver
module load netcdf/$netcdf_ver

module list

mkdir -p ../exec

cd ./nsst.fd
./makefile.sh
