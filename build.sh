#! /usr/bin/env bash
#
# Author: Larissa Reames CIWRO/NOAA/NSSL/FRDD

set -eux

target=${1:-"NULL"}
compiler=${compiler:-"intel"}
echo $target, $compiler
if [[ "$target" == "linux.*" || "$target" == "macosx.*" ]]; then
    unset -f module
    set +x
    source ./modulefiles/build.$target > /dev/null
    set -x
else
    set +x
   #source ./machine-setup.sh
   #module use ./modulefiles
   #module load build.$target.$compiler.lua
   #module list
    set -x
fi

target="derecho"

if [[ "$target" == "cheyenne" ]] ; then # $target set in ./machine-setup.sh
   module purge
   module restore modules_intel19 # from user schwartz
fi
if [[ "$target" == "derecho" ]] ; then
   module --force purge
   module restore default # default is from user schwartz
   module swap hdf5/1.12.2 hdf5-mpi/1.12.2
   module swap netcdf      netcdf-mpi
   module load parallel-netcdf/1.12.3
   module load parallelio/2.6.2
   module list
fi

if [[ "$target" == "hera" || "$target" == "orion" || "$target" == "wcoss2" ]]; then
  #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Debug"
   CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF"
   #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DENABLE_DOCS=ON -DBUILD_TESTING=ON"
else
  #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Debug"
  CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF"
fi

rm -fr ./build
mkdir ./build && cd ./build

cmake .. ${CMAKE_FLAGS}

make -j 8 VERBOSE=1
make install

#make test
#ctest -I 4,5

exit
