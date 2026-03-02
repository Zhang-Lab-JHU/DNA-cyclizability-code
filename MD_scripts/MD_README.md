vmd_script.txt contains instructions for generating input files for MD simulation in NAMD from nucleotide sequence

most_cyclizable_initial_files folder contains pre-generated files for NAMD simulations for the most cyclizable sequence constructed from the correlation function model

MD pipeline consists of:
  1. Minimization - min_tut.conf
  2. Constant-volume equilibration - equil1.conf
  3. Constant-pressure equilibration - equil2.conf
  4. Final simulation run - product.conf
which should be run sequentially

par_all36_na.prm and toppar_water_ions.mod.str are parameter files for NAMD
