for nx in 64 128 256; do
	echo $nx
	sed -i "/nx = /c\   integer,  parameter :: nx = $nx      ! Number of grid points (excluding boundaries)." example/exptA_GL/ginzburg_landau_base.f90
	grep nx example/exptA_GL/ginzburg_landau_base.f90
	fpm run --example main
done
