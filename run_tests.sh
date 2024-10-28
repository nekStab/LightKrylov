# Gfortran.
python deployment.py

fpm test --compiler "gfortran" --flag "-O0 -g3 -fbacktrace -Wall -Wextra -fcheck=all -pedantic -Wconversion -fbounds-check -ffpe-trap=zero,overflow" --verbose
fpm test --compiler "gfortran" --profile release --flag "-march=native -mtune=native"

# Ifort
fpm test --compiler "ifort" --flag "-O0 -g -traceback -check all -check bounds -check all -check uninit, -ftrapuv -debug all" --verbose
fpm test --compiler "ifort" --profile release --flag "-xHost"

# Ifx
#fpm test --compiler "ifx"
#fpm test --compiler "ifx" --profile release --flag "-O3 -xHost"
