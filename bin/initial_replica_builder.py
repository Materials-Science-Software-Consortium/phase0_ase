#!/usr/bin/env python

import optparse
import sys

from ase import io
from ase.neb import NEB

def run():
    parser = optparse.OptionParser(usage='%prog [options]',description='create initial replicas for the PHASE/0 NEB feature')
    parser.add_option('-i','--initial',dest='initial_file',type=str,default=None,help='initial coordinates. Required')
    parser.add_option('-f','--final',dest='final_file',type=str,default=None,help='final coordinates. Required')
    parser.add_option('--initial_format',dest='initial_format',type=str,default='phase0-out',help='format for the initial coordinate file. Defaults to phase0-out')
    parser.add_option('--final_format',dest='final_format',type=str,default='phase0-out',help='format for the final coordinate file. Defaults to phase0-out')
    parser.add_option('-n','--nimage',dest='nimage',type=int,default=6,help='number of intermediate images. defaults to 6')
    parser.add_option('-l','--linear',dest='linear',action='store_true',default=False,help='specify this option in order to create replicas by simple linear interpolation. the default is to create replicas by the IDPP method.')

    (options,args) = parser.parse_args()

    if options.initial_file is None or options.final_file is None:
        print '--initial and --final are required options'
        sys.exit()
    initial = io.read(options.initial_file,format=options.initial_format,wrap=False)
    final = io.read(options.final_file,format=options.final_format,wrap=False)

    images = []
    images.append(initial)
    for i in range(options.nimage):
        images.append(initial.copy())
    images.append(final)

    neb = NEB(images)
    neb.interpolate()
    if not options.linear:
        neb.idpp_interpolate()

    for i in range(options.nimage+2):
        image = neb.images[i]
        natm = len(image)
        im = open('image'+str(i)+'.data','w')
        im.write('coordinate_system=cartesian\n')
        im.write('\n')
        im.write('#units angstrom\n')
        im.write('\n')
        symb = image.get_chemical_symbols()
        spos = image.get_positions()
        for iat in range(natm):
            im.write(str(symb[iat])+' '+str(spos[iat][0])+' '+str(spos[iat][1])+' '+str(spos[iat][2])+'\n')
        im.close()
    io.write('initial_replicas.xyz',images) 
    
if __name__ == '__main__':
    run()

