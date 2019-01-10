from ase.db import connect
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
from ase.visualize import view
from ase.io.trajectory import Trajectory
from gpaw import GPAW, PW, FermiDirac, restart
import sys
import os
from ase.io import read, write

traj = Trajectory('fesi_8atoms_id10.traj')
atoms = traj[-1]
write('fesi_system.xyz', atoms)
