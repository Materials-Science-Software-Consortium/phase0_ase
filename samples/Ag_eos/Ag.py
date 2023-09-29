#!/usr/bin/env python

from ase import Atoms
from ase.calculators.phase0 import phase0
from ase.units import kJ
from ase.eos import calculate_eos

a = 4.15  # approximate lattice constant
b = a / 2
ag = Atoms('Ag',
           cell=[(0, b, b), (b, 0, b), (b, b, 0)],
           pbc=1,
           calculator=phase0(ecut=25.0))
eos = calculate_eos(ag,10,0.1,trajectory='Ag.traj')
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
eos.plot('Ag-eos.png')

