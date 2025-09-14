#for nx in 64 128 256 512; do
for nx in 512 1024; do
	echo $nx
	sed -i "/nx = /c\   integer,  parameter :: nx = $nx      ! Number of grid points (excluding boundaries)." example/kexptA_GL/ginzburg_landau_base.f90
	fpm run --example main | tee kexpm_test.txt
	mv kexpm_test.txt kexpm_test_nx_$nx.txt
done
