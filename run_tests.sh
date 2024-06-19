# Gfortran.
python deployment.py

fpm test --compiler "gfortran" --flag "-O0 -g3 -fbacktrace -Wall -Wextra -fcheck=all -pedantic -Wconversion -fbounds-check -ffpe-trap=zero,overflow" --verbose
fpm test --compiler "gfortran" --flag "-O3 -march=native -mtune=native"

# Ifort
#fpm test --compiler "ifort"
#fpm test --compiler "ifort" --flag "-O3 -xhost"

# Ifx
#fpm test --compiler "ifx"
#fpm test --compiler "ifx" --flag "-O3 -xhost"
