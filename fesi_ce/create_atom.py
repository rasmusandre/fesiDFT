from ase.build import bulk
from ase.visualize import view

build_atoms = bulk('Fe', 'bcc', a=2.876, cubic=True)*(2,2,2)
view(build_atoms)
