#!/usr/bin/env python

import optparse
import os
import sys

from ase import Atoms
from ase.io import read,write
from ase.calculators.phase0 import phase0
from ase.optimize import QuasiNewton
from ase.eos import calculate_eos
from ase.units import kJ

def run():
    
    parser = optparse.OptionParser(usage="%prog [options]",description="PHASE/0 wrapper script.\
 environment variable ASE_PHASE0_COMMAND and PHASE_PP_PATH must be set prior to the execution of this script.\
 ASE_PHASE0_COMMAND specifies the command used to boot PHASE/0.\
 example : 'mpirun -np 8 phase ne=2 nk=4'.\
 PHASE_PP_PATH specifies the directory under which the PHASE/0 pseudopotential files reside.\
 example : $HOME/pp ")
    parser.add_option('-f','--file',dest='coord_file',type=str,default=None,help='the file name for the coordinate file. Required.')
    parser.add_option('-i','--index',dest='index',type=int,default=-1, 
                      help='the index of the coordinate if the specified file supports frame. defaults to -1, i.e., the last frame ')
    parser.add_option('-t','--type',dest='coord_type',type=str,default=None,
                      help='type of the coordinate file. defaults to None, in which case the ase resolver resolves the file type' )
    parser.add_option('-e','--ecut',dest='ecut',type=float,default=25.0,
                      help='cutoff energy for the WFs in Rydberg units. the default value is 25')
    parser.add_option('-k','--kpts',dest='kpts',type=str,default='3.5',
                      help='specify the kpoint sampling. specify either one float value or three comma-seperatred integers. \
                            the former specifies the density of the kpoints in A^-1 units.\
                            the latter corresponds to the mesh for the a,b and c axis. the default value is 3.5')
    parser.add_option('-x','--fmax',dest='fmax',type=float,default=0.01,help='convergence criterion (max. force) in eV/A units. \
                       the default value is 0.01')
    parser.add_option('-s','--spin',dest='spin',action='store_true',default=False, 
                      help='enable this option in order to perform spin-polarized calculations.')
    parser.add_option('','--nosymm',dest='nosymm',action='store_true',default=False, 
                      help='enable this option in order to disable symmetry.')
    parser.add_option('-y','--symmetry_check_criterion',dest='symmetry_check_criterion',type=float,default=None, 
                      help='specify the criterion for symmetry check. defaults to 1e-12')
    parser.add_option('','--inversion',dest='inversion',action='store_true',default=False, 
                      help='enable this option in order to set the sw_inversion parameter to \'on\'.')
    parser.add_option('-z','--zeta',dest='zeta',type=str,default=None,help='specify the initial per-element spin polarization.\
                      use the following format : --zeta=elem1:zeta1,elem2:zeta2 ...')
    parser.add_option('-p','--pp',dest='pp',default=None,\
                       help='specify pseudopotential files in the following format : --pp=elem1:pp1,elem2:pp2,...')
    parser.add_option('-c','--clean',dest='clean',action='store_false',default=False)
    parser.add_option('-m','--mode',dest='mode',type="choice",choices=['optimize','eos','preproc','input_only'],default='optimize',\
                       help='specify the operation mode. one of : optimize, eos, preproc, or input_only. \
 optimize : structural optimization by the ase quasi newton solver (default) eos      : equation of state calculation\
 preproc  : run PHASE/0 in preparation mode\
 input_only : create input and exit')
    parser.add_option('-d','--delta',dest='delta',type=float,default=0.1,
                       help='the eps value for the calculate_eos method. does not make sense when --mode is not eos')
    parser.add_option('-n','--npoints',dest='npoints',type=int,default=10,
                       help='the number of points for the calculate_eos method. does not make sense when --mode is not eos')

    (options,args) = parser.parse_args()

    comm = os.getenv('ASE_PHASE0_COMMAND')
    if comm is None:
        print('environment variable ASE_PHASE0_COMMAND is unset')
        sys.exit()
    pp = os.getenv('PHASE_PP_PATH')
    if pp is None:
        print('environment variable PHASE_PP_PATH is unset')
        sys.exit()

    if options.coord_file is None:
        print('coordinate file not specified')
        sys.exit()

    pref = options.coord_file.split('.')[0]

    atoms = read(options.coord_file,index=options.index,format=options.coord_type)
    kpts = None
    kps = options.kpts.strip().split(',')
    if len(kps)<3:
        kpts = float(kps[0])
    else:
        kpts = (int(kps[0]),int(kps[1]),int(kps[2]))

    zeta = None
    if options.zeta is not None:
        z = options.zeta.split(',')
        zeta = {}
        for zz in z:
            zzz = zz.split(':')
            zeta[zzz[0]] = zzz[1]
    pp = None
    if options.pp is not None:
        p = options.pp.split(',')
        pp = {}
        for pppp in p:
            ppp = pppp.split(':')
            pp[ppp[0]] = ppp[1]
    symm = not options.nosymm
    phase = phase0(ecut=options.ecut,kpts=kpts,spin=options.spin,zeta=zeta,pp=pp, \
                   symmetry=symm,inversion=options.inversion,preproc=options.mode == 'preproc', \
                   symmetry_check_criterion=options.symmetry_check_criterion)
    atoms.set_calculator(phase)

    if options.mode == 'optimize':
        qn = QuasiNewton(atoms, trajectory=pref+'.traj')
        qn.run(fmax=options.fmax)
    elif options.mode == 'eos':
        eos = calculate_eos(atoms,npoints=options.npoints,eps=options.delta,trajectory=pref+'.traj')
        v,e,B = eos.fit()
        print('equilibrium volume : '+str(v) +' A^3')
        print('equilibrium energy : '+str(e) + ' eV')
        print('bulk modulous : '+str(B / kJ * 1.0e24) + ' GPa')
        eos.plot(pref+'.png')
    elif options.mode == 'preproc':
        phase.export_input(atoms)
        os.system(os.getenv('ASE_PHASE0_COMMAND'))
    elif options.mode == 'input_only':
        phase.export_input(atoms)

    if options.clean:
        phase.clean()
        
if __name__ == '__main__':
    run()

