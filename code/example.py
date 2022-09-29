import os
import sys
import pathlib as pl
import shutil

sys.path.append(os.path.abspath(r'../../'))
from utilities import PointCloudDataGen, classes

if __name__ == '__main__':
	cwd = pl.Path(os.getcwd())

	# ASC surfaces
	asc1 = cwd / 'sample1.asc'
	asc2 = cwd / 'sample2.asc'
	PointCloudDataGen.spherical_asperities(6.67*2, 25, 6.67, asc1)
	PointCloudDataGen.planar(60, 240, asc2)

	# create simulation directory
	simulation_dir = cwd / "simulation"
	if simulation_dir.exists():
		shutil.rmtree(simulation_dir)
	os.mkdir(simulation_dir)

	# surface input files
	write_file1 = simulation_dir / "sample1.inp"
	write_file2 = simulation_dir / "sample2.inp"

	surface1 = classes.ASC_surface(asc1)
	surface2 = classes.ASC_surface(asc2)
	# classes.surfaces_reinterpolation(surface1, surface2)

	volume1 = classes.VolumeMesh(surface1)
	volume1.write_inp_file(write_file1)

	volume2 = classes.VolumeMesh(surface2, layers=20)
	volume2.write_inp_file(write_file2)
