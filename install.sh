#!/bin/bash

FC="mpif90"
FFLAGS=""

if command -v mpiifort >/dev/null 2>&1; then
    echo "Compiler found: $(mpiifort -v)"
    FC="mpiifort"
    FFLAGS+=" -O3 -no-prec-div -fp-model fast=2 -xHost -g -traceback" 
   if command -v mpiifx >/dev/null 2>&1; then
        echo "Using LLVM version of Intel compiler..."
        FC="mpiifx"
        FFLAGS+=" -qmkl"
   else
        echo "LLVM not found, using classic."
        FFLAGS+=" -ipo -qmkl"
   fi
else
    echo "Compiler found: $(mpif90 --version | head -n1)"
    FFLAGS+=" -O3 -march=native -funroll-loops -ffast-math -g -fbacktrace"
fi
#FFLAGS+=" -fanalyser" # Enable the Clang Static Analyzer for code analysis
# export FFLAGS
# export FC

echo "Compiler: $FC"
echo "Flags: $FFLAGS"

if ! command -v fpm &> /dev/null
then
   echo "fpm could not be found"
   exit
fi

fpm clean
read -p "Do you want to run tests before installing? (y/n) " answer
case ${answer:0:1} in
   y|Y )
      fpm test --verbose --compiler "$FC" --flag "$FFLAGS" || { echo 'fpm test failed' ; exit 1; }
   ;;
   * )
      echo "Skipping tests..."
   ;;
esac
fpm install --verbose --compiler="$FC" --flag="$FFLAGS" || { echo 'fpm install failed' ; exit 1; }


# fpm test --compiler 'mpif90' --flag "-O3 -march=native -funroll-loops -ffast-math -g -fbacktrace"
# fpm test --compiler 'mpiifort' --flag "-O3 -no-prec-div -fp-model fast=2 -xHost -g -traceback -ipo -mkl"
# fpm test --compiler 'mpiifx' --flag " -O3 -no-prec-div -fp-model fast=2 -xHost -g -traceback -qmkl"