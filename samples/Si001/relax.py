#!/usr/bin/env python

from ase.calculators.phase0 import phase0
from ase.optimize import QuasiNewton
from ase.io import read,write

atoms = read('Si001.cif')
phase = phase0(ecut=16,kpts=(4,8,1))
atoms.set_calculator(phase)

q = QuasiNewton(atoms,trajectory='Si001.traj')
q.run(fmax=0.005)

traj = read('Si001.traj',index=':') # read all frames by specifying index=':'
write('Si001_relax.cif',traj)

