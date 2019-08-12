all:
	python -m numpy.f2py -c weno.f90 -m weno
