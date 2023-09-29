"""
This module contains functionality for reading and writing an ASE
Atoms object in PHASE/0 nfinp.data and nfdynm.data format

"""
from ase.utils import basestring
from ase.data import atomic_masses,chemical_symbols,atomic_numbers
from ase import Atoms, Atom

import sys
import os
from os.path import join,exists

import re

import logging
logformat  = '%(levelname)10s: %(message)s'
logtimefmt = '%a %d %b %H:%M:%S'
logging.basicConfig(level=logging.INFO,
                    format=logformat,
                    datefmt=logtimefmt
                    )
logger = logging.getLogger('phase0')

Hartree_2_eV = 27.21139615
Bohr_2_Angstrom = 0.5291772480
Hartree_Bohr_2_eV_Angstrom = Hartree_2_eV/Bohr_2_Angstrom
Angstrom_2_Bohr = 1.0/Bohr_2_Angstrom
Hartree_Bohr3_2_eV_Ang3 = Hartree_2_eV/Bohr_2_Angstrom/Bohr_2_Angstrom/Bohr_2_Angstrom

BLOCK = 'block'
TABLE = 'table'
PRIMITIVE = 'primitive'
VECTOR = 'vector'

FRACTIONAL = 0
CARTESIAN = 1
BOHR = 0
ANGSTROM = 1
NM = 2

unit_type            = 'unit_type'
unit_type_none       = 'none'
unit_type_length     = 'length'
unit_type_energy     = 'energy'
unit_type_time       = 'time'
unit_type_longtime   = 'longtime'
unit_type_velocity   = 'velocity'
unit_type_force      = 'force'
unit_type_pressure   = 'pressure'
unit_type_mass       = 'mass'
unit_type_angle      = 'angle'

unit_spec={
unit_type_none:[''],
unit_type_length:['bohr','angstrom','nm'],
unit_type_energy:['hartree','ev','rydberg'],
unit_type_time:['au_time','fs','ps','ns'],
unit_type_longtime:['sec','min','hour','day'],
unit_type_velocity:['bohr/au_time','bohr/fs','angstrom/fs','angstrom/au_time','nm/fs','nm/au_time'],
unit_type_force:['hartree/bohr','hartree/angstrom','hartree/nm','ev/angstrom','ev/bohr','ev/nm','rydberg/bohr','rydberg/angstrom','rydberg/nm'],
unit_type_pressure:['hartree/bohr3','hartree/angstrom3','hartree/nm3','ev/angstrom3','ev/bohr3','ev/nm3','rydberg/angstrom3','rydberg/bohr3','rydberg/nm3'],
unit_type_mass:['au_mass','atomic_mass'],
unit_type_angle:['degree','radian']
}

output000_p = re.compile('output\d\d\d')
jobstatus000_p = re.compile('jobstatus\d\d\d')
def restartable(directory='.'):
    contfile_list = ['zaj.data','nfchgt.data','continue.data','continue_bin.data']
    files = os.listdir(directory)
    for cfile in contfile_list:
        exi = False
        for f in files:
            if cfile == f:
                exi = True
                continue
        if not exi:
            return False
    return True
    
def get_outputxxx(directory='.'):
    files = os.listdir(directory)
    ret = []
    for f in files:
        if output000_p.match(f) is not None:
            ret.append(f)
    return ret

def get_jobstatusxxx(directory='.'):
    files = os.listdir(directory)
    ret = []
    for f in files:
        if jobstatus000_p.match(f) is not None:
            ret.append(f)
    return ret

def get_mostrescent_outputxxx(directory='.'):
    """ get the most rescent PHASE/0 log file (outputxxx file) under the specified directory"""
    files = os.listdir(directory)
    ret = None
    mod = -1
    for f in files:
        if output000_p.match(f) is not None:
            t = os.path.getmtime(f)
            if t>mod:
                ret = f
    return ret
    
def get_default_ppfile_from_dir(element_name,dire=None,relative=False,full_path=False):
    ''' obtain the path to the default pp file corresponding to the specified element. 
        will return None if the pp file could not be resolved. if relative is set to True,
        will return the file name by relative path with respect to cwd. '''
    if dire is None:
        if 'PHASE_PP_PATH' in os.environ:
            dire = os.environ['PHASE_PP_PATH']
        else:
            logger.error("environment variable PHASE_PP_PATH unset")
            sys.exit()
    if dire is None:
        return None
    dirlist = os.listdir(dire)
    pps=[]
    if element_name is None:
        return None
    for file in dirlist:
        if file.startswith(element_name+"_"):
            ppfile_elem = file.split('.')[0].split('_')
            id = ppfile_elem[len(ppfile_elem)-1]
            gga_lda=file.split('_ggapbe_')
            if len(gga_lda)>=2 and file.endswith('.pp') and len(id)<=2:
                pps.append(file)
    if len(pps)>=1:
        if not full_path:
            return pps[len(pps)-1]
        else:
            return dire + '/' + pps[len(pps)-1]
    return None

def get_pp_path_as_list(input):
    """ obtain a list of pp files. Will return None if atomic coordinates 
        and element info are  undefined in the input """
    atom_table = input.get_root_entry().get_entry('structure.atom_list.atoms.table')
    elem_table = input.get_root_entry().get_entry('structure.element_list.table')

    if atom_table is None:
        logger.error('atom table is undefined.')
        return None

    nat = atom_table.get_num_data()
    elem_dict={}
    for i in range(nat):
        elem = atom_table.get_data(i,'element')
        if not elem in elem_dict.keys():
            elem_dict[elem] = False

    nelem = elem_table.get_num_data()
    ppfiles = []
    for i in range(nelem):
        elemname = elem_table.get_data(i,'element')
        if elemname is None:
            logger.error("element name undefined")
            return None
        belename = re.sub("\d","",elemname)
        ppfiles.append(get_default_ppfile_from_dir(belename))
    return ppfiles

def get_ws(length):
    ''' return white-space characters '''
    ret =''
    for i in range(length):
        if i==0:
            continue
        ret += '  '
    return ret

def get_value(spec,val,index=0):
    ''' obtain converted value according to spec.
        if spec ==None, or if any errors are encountered,
        will return str(val). '''
    if spec==None or not val_type in spec:
        return str(val)

    vtyp = spec[val_type].split(',')[index].strip()
    if vtyp == val_type_float:
        try:
            return float(val)
        except ValueError:
            return str(val)
    if vtyp == val_type_int:
        try:
            return int(val)
        except ValueError:
            return str(val)
    if vtyp == val_type_bool:
        return pyutil.parse_bool(val)
    return str(val)

def resolve_unit(unit_name):
    ''' return the 'type' of the unit associated with the given unit name.
        for example, 'bohr' will return 'length'. will return None if
        no such unit exists.
        parameter : unit_name the name of the unit. '''
    if unit_name == None:
        return None
    for typ,units in unit_spec.items():
        if unit_name.lower() in units:
            return typ
    return None

def read_phase0(filename='nfinp.data',wrap=False):
    if isinstance(filename, basestring):
        f = open(filename)
    else:  # Assume it's a file-like object
        f = filename
    fname = filename.name
    phasei = Input(coord_file=fname)
    return phasei.get_atomic_configuration(wrap)

def write_phase0(filename, atoms, cartesian=False, angstrom=False, mobile=False):
    if isinstance(filename, basestring):
        f = open(filename, 'w')
    else:  # Assume it's a 'file-like object'
        f = filename.name
    if isinstance(atoms, (list, tuple)):
        if len(atoms) > 1:
            raise RuntimeError('Don\'t know how to save more than ' +
                               'one image to PHASE/0 input')
        else:
            atoms = atoms[0]
    phasei = Input(coord_file=f)
    phasei.swap_atomic_configuration(atoms,to_file=f, cartesian=cartesian, angstrom=angstrom,mobile=mobile)
    phasei.save()

def update_continue_data(atoms,directory='.',fname='continue.data'):
    fcoords = atoms.get_scaled_positions()
    ccoords = atoms.get_positions() 
    cont = open(directory+'/'+fname)
    newcont = open(directory+'/conttmp','w')
    natm = -1 
    in_isol = False
    in_natm = False
    in_pos = False
    in_cps = False
    in_conver = False
    in_iter = False
    atm_count = 0
    for line in cont:
        sline = line.strip()

        if in_iter:
            words = sline.split()
            words[2] = '0'
            newcont.write(words[0]+'  '+words[1]+'  '+words[2]+'\n')
            in_iter = False
            continue
        if in_pos:
            coor = fcoords[atm_count]
            newcont.write(str(coor[0])+' '+str(coor[1])+' '+str(coor[2])+'\n')
            atm_count += 1
            if atm_count == natm:
                in_pos = False
            continue

        if in_cps:
            coor = ccoords[atm_count]
            newcont.write(str(coor[0]*Angstrom_2_Bohr)+' '+str(coor[1]*Angstrom_2_Bohr)
                     +' '+str(coor[2]*Angstrom_2_Bohr)+'\n')
            atm_count += 1
            if atm_count == natm:
                in_cps = False
            continue

        if in_isol:
            newcont.write('         56\n') # Kosugi solver
            in_isol = False
            continue
        if in_conver:
            newcont.write('         0\n')
            in_conver = False
            continue

        newcont.write(line)

        if sline.startswith('iteration'):
            in_iter = True
            continue

        if sline.startswith('(natm)'):
            in_natm = True
            continue

        if in_natm:
            natm = int(line)
            in_natm = False
            continue

        if sline.startswith('(pos)'):
            in_pos = True
            atm_count = 0
            continue

        if sline.startswith('(cps)'):
            in_cps = True
            atm_count = 0
            continue

        if sline.startswith('convergence'):
            in_conver = True
            continue
        if sline.startswith('isolver'):
            in_isol = True
            continue
    newcont.close()
    cont.close() 
    os.remove(fname)
    os.rename('conttmp','continue.data')

def read_phase0_out(filename='nfdynm.data', index=-1, wrap=True, fname_from_fnamesdata=False):
    """ Import PHASE/0 nfdynm.data type file.
    Reads unitcell, atom positions, and forces from the nfdynm.data 
    (file specified by the F_DYNM pointer) file
    if present, reads the energies from the nfefn.data
    (file specified by the F_ENF pointer) file
    """
    from ase.calculators.singlepoint import SinglePointCalculator
    import numpy as np

    nfefndata = 'nfefn.data'
    if fname_from_fnamesdata:
        fnames = FileNamesData()
        filename = fnames.get_file_name('F_DYNM')
        nfefndata = fnames.get_file_name('F_ENF')
    energies = []
    if exists(nfefndata):
        nfefn = open(nfefndata)
        for line in nfefn:
            line = line.strip()
            if line.startswith('iter') or len(line)==0:
                continue
            words = line.split()
            if len(words) < 3:
                print('invalid line found in the F_ENF file')
            try:
                energies.append(float(words[2])*Hartree_2_eV)
            except ValueError:
                pass
        nfefn.close()

    IN_HEADER=0
    IN_DATA=2

    if isinstance(filename, basestring):
        f = open(filename)
    else:  # Assume it's a file-like object
        f = filename
    data = f.readlines()
    ntyp=0
    natm=0
    avec=[]
    bvec=[]
    cvec=[]
    atm_typ=[]
    species_map = {}
    symbol = []
    mode = IN_HEADER
    atm_count=0
    images = []
    forcess = []
    atoms = None
    cell = None
    for lineno,line in enumerate(data):
        line = line.strip()
        if mode==IN_HEADER and line.startswith('cps'):
            for at in atm_typ:
                symbol.append(species_map[at])
            cell = (avec,bvec,cvec)
            mode = IN_DATA
            atm_count=0

        if line.startswith('cps'):
            atoms = Atoms(pbc=wrap)
            forces = []
            forcess.append(forces)
            atoms.set_cell(cell)
            images += [atoms]
            atm_count=0
            continue

        if line.startswith('#'):
            mode=IN_HEADER
            words = line.split('#')
            if len(words)<2:
                continue
            words = line.split('#')[1].split()
            if len(words)==0:
                continue
            word0 = words[0]
            try:
                if word0=='a_vector':
                    avec=[]
                    bvec=[]
                    cvec=[]
                    atm_typ=[]
                    species_map={}
                    symbol=[]
                    avec.append(float(words[2])*Bohr_2_Angstrom)
                    avec.append(float(words[3])*Bohr_2_Angstrom)
                    avec.append(float(words[4])*Bohr_2_Angstrom)
                elif word0=='b_vector':
                    bvec.append(float(words[2])*Bohr_2_Angstrom)
                    bvec.append(float(words[3])*Bohr_2_Angstrom)
                    bvec.append(float(words[4])*Bohr_2_Angstrom)
                elif word0=='c_vector':
                    cvec.append(float(words[2])*Bohr_2_Angstrom)
                    cvec.append(float(words[3])*Bohr_2_Angstrom)
                    cvec.append(float(words[4])*Bohr_2_Angstrom)
                elif word0=='ntyp':
                    ntyp = int(words[2])
                    natm = int(words[5])
                elif word0=='(natm->type)':
                    for word in words:
                        if word==word0:
                            continue
                        atm_typ.append(int(word))
                elif word0.startswith('(speciesname)'):
                    species_map[int(words[1])] = words[3]
            except ValueError:
                raise ValueError("encountered error at line : "+str(lineno)+" in file "+filename.name)

        if mode==IN_DATA:
            words = line.split()
            try:
                atoms += Atom(species_map[atm_typ[atm_count]],
                [float(words[1])*Bohr_2_Angstrom,float(words[2])*Bohr_2_Angstrom,float(words[3])*Bohr_2_Angstrom])
                forces.append(np.array([float(words[4])*Hartree_Bohr_2_eV_Angstrom,
                           float(words[5])*Hartree_Bohr_2_eV_Angstrom,
                           float(words[6])*Hartree_Bohr_2_eV_Angstrom]))
                atm_count+=1
            except ValueError:
                raise ValueError("encountered error at line : "+str(lineno)+" in file "+filename.name)
    for i,atoms in enumerate(images):
        ene = 0.0
        forces = []
        for ii in range(len(atoms)):
            forces.append([0.0,0.0,0.0])
        if i<len(energies):
            ene = energies[i]
        if i<len(forcess):
            forces = forcess[i]
        atoms.set_calculator(SinglePointCalculator(atoms,energy=ene,forces=forces))

    # return requested images, code borrowed from ase/io/vasp.py
    if isinstance(index, int):
        return images[index]
    else:
        step = index.step or 1
        if step > 0:
            start = index.start or 0
            if start < 0:
                start += len(images)
            stop = index.stop or len(images)
            if stop < 0:
                stop += len(images)
        else:
            if index.start is None:
                start = len(images) - 1
            else:
                start = index.start
                if start < 0:
                    start += len(images)
            if index.stop is None:
                stop = -1
            else:
                stop = index.stop
                if stop < 0:
                    stop += len(images)
        return [images[i] for i in range(start, stop, step)]

class Line:
    ''' the class which encapsulates a line of an input file'''
    def __init__(self,line,line_no):
        self._line=line
        self._line_no=line_no

    def get_line_no(self):
        ''' get the no. of the line '''
        return self._line_no

    def get_line(self):
        ''' get the line it self '''
        return self._line

    def set_line(self,line):
        ''' the setter for the line '''
        self._line=line

    def add_line(self,line):
        ''' the line may also be 'appended' '''
        self._line += line
    
    def __str__(self):
        return str(self._line_no)+":"+self._line

class Input():
    ''' class which represents a PHASE/0 input file

        typical usage are as follows...

        inp = input.Input()

        root_block=inp.get_root_entry() -> obtain the 'root entry' of F_INP.
        control = root_block.get_entry('control') -> obtain the 'control' block
        control.is_block() -> True
        cpuma=control.get_entry('cpumax') -> get an instance of class PrimitiveEntry
        print cpuma.get_value()+cpuma.get_unit() -> get something like 1 day
        cpuma.set_value('2') -> the 'control.cpumax' parameter will become '2'
        ...
        ...
        ...
        input.save() -> saves to F_INP
        input.save_to(file) -> saves to file '''
    root_block=None
    curr_pblock=None
    ljust_len=70
    def __init__(self,specfile_name=None,specdir=None,coord_file='nfinp.data'):
        ''' parameter : specfile_name the filename of the specification file.
            if given, type conversions of primitive types will be performed.
            parameter : specdir the directory on which the spec file resides.
            defaults to psf.data.inputspec'''
        self.nwarn=0
        self.nerror=0
        self.inp_exists=False
        self.filename_data = FileNamesData(self,nfinp=coord_file)
        if coord_file is None:
            if not self.filename_data.file_names_data_exists():
                return    

        self.filename_data.debug_output()
        if coord_file is None:
            self.file_name = self.filename_data.get_file_name('F_INP')
        else:
            self.file_name = coord_file
        if self.file_name is None:
            return

        if not os.path.exists(self.file_name):
            self.inp_exists = False
        else:
            self.inp_exists = True

        self.root_block = BlockEntry("root",self,0)
        self.curr_pblock = self.root_block
        if specfile_name!=None:
            self.__init_spec(specfile_name,specdir)

        self.__parse_input()
        self.__dump()

    def get_atomic_configuration(self,wrap=False):
        if (self.get_root_entry().get_entry('structure.atom_list') is None or 
           self.get_root_entry().get_entry('structure.atom_list.atoms') is None):
            logger.error('structure.atom_list.atoms undefined. ')
            return None

        coordsys = self.get_root_entry().get_entry('structure.atom_list.coordinate_system')
        cs = FRACTIONAL
        if coordsys is not None and coordsys.get_value()=='cartesian':
            cs = CARTESIAN

        uni_coords = Bohr_2_Angstrom
        if self.get_root_entry().get_entry('structure.atom_list.atoms').get_default_unit('length')\
           is not None:
            struni = self.get_root_entry().get_entry('structure.atom_list.atoms').get_default_unit('length')
            if struni.lower() == 'angstrom':
                uni_coords = 1.0
            elif struni.lower()=='nm':
                uni_coords = 10.0
        
        cell = None
        a_vector = self.get_root_entry().get_entry('structure.unit_cell.a_vector')
        b_vector = self.get_root_entry().get_entry('structure.unit_cell.b_vector')
        c_vector = self.get_root_entry().get_entry('structure.unit_cell.c_vector')
        celluni = Bohr_2_Angstrom
        cellu = self.get_root_entry().get_entry('structure.unit_cell').get_default_unit('length')
        if cellu is not None and cellu.lower()!='bohr':
            if cellu.lower() == 'angstrom':
                celluni = 1.0
            elif cellu.lower()=='nm':
                celluni = 10.0
        if a_vector is not None and b_vector is not None and c_vector is not None:
            a0 = float(a_vector.get_value()[0])*celluni
            a1 = float(a_vector.get_value()[1])*celluni
            a2 = float(a_vector.get_value()[2])*celluni
            b0 = float(b_vector.get_value()[0])*celluni
            b1 = float(b_vector.get_value()[1])*celluni
            b2 = float(b_vector.get_value()[2])*celluni
            c0 = float(c_vector.get_value()[0])*celluni
            c1 = float(c_vector.get_value()[1])*celluni
            c2 = float(c_vector.get_value()[2])*celluni
            cell = ((a0,a1,a2),(b0,b1,b2),(c0,c1,c2))

        if cell is None:
            a = float(self.get_root_entry().get_entry('structure.unit_cell.a').get_value())*celluni
            b = float(self.get_root_entry().get_entry('structure.unit_cell.b').get_value())*celluni
            c = float(self.get_root_entry().get_entry('structure.unit_cell.c').get_value())*celluni
            alpha = float(self.get_root_entry().get_entry('structure.unit_cell.alpha').get_value())*celluni
            beta  = float(self.get_root_entry().get_entry('structure.unit_cell.beta').get_value())*celluni
            gamma = float(self.get_root_entry().get_entry('structure.unit_cell.gamma').get_value())*celluni
            if a is not None and b is not None and c is not None and \
               alpha is not None and beta is not None and gamma is not None:
               cell = (a,b,c,alpha,beta,gamma)

        atoms = Atoms(pbc=wrap)
        atoms.set_cell(cell)
        (avec,bvec,cvec) = atoms.get_cell() 
        atm_table = self.get_root_entry().get_entry('structure.atom_list.atoms.table')
        natm = atm_table.get_num_data()
        for i in range(natm):
            elem = atm_table.get_data(i,'element')
            rx = float(atm_table.get_data(i,'rx'))
            ry = float(atm_table.get_data(i,'ry'))
            rz = float(atm_table.get_data(i,'rz'))
            xx = 0.0
            yy = 0.0
            zz = 0.0
            if cs == FRACTIONAL:
               xx = avec[0]*rx+bvec[0]*ry+cvec[0]*rz
               yy = avec[1]*rx+bvec[1]*ry+cvec[1]*rz
               zz = avec[2]*rx+bvec[2]*ry+cvec[2]*rz
            elif cs == CARTESIAN:
               xx = rx * uni_coords
               yy = ry * uni_coords
               zz = rz * uni_coords
            atoms += Atom(elem,[xx,yy,zz])

        return atoms

    structure_template=\
    "\
structure{\n\
    atom_list{\n\
        atoms{\n\
            #tag element rx ry rz\n\
        }\n\
    }\n\
    unit_cell{\n\
    }\n\
    element_list{\n\
        #tag element atomicnumber mass\n\
    }\n\
}\n\
    "
    def swap_atomic_configuration(self,atoms,to_file=None, cartesian=False, angstrom=False, mobile=False,zeta=None):
        if not cartesian:
            coord = atoms.get_scaled_positions()
        else:
            coord = atoms.get_positions()
        coord_fac = 1.0
        if not angstrom and cartesian:
            coord_fac = Angstrom_2_Bohr
        symbols = atoms.get_chemical_symbols()
        natms = len(symbols)
        (avec,bvec,cvec) = atoms.get_cell() 
        cell_fac = 1.0
        if not angstrom:
            cell_fac = Angstrom_2_Bohr 

        atomtable = self.get_root_entry().get_entry('structure.atom_list.atoms.table')
        blunitcell = self.get_root_entry().get_entry('structure.unit_cell')
        elem_table = self.get_root_entry().get_entry('structure.element_list.table')
        if atomtable is None or blunitcell is None or elem_table is None and to_file is not None:
            file=open(to_file,'w')
            file.write(self.structure_template)
            file.flush()
            file.close()
            self.__parse_input(to_file)
            atomtable = self.get_root_entry().get_entry('structure.atom_list.atoms.table')
            blunitcell = self.get_root_entry().get_entry('structure.unit_cell')
            elem_table = self.get_root_entry().get_entry('structure.element_list.table')
        atomtable.reset_data()
        elem_table.reset_data()
        self.get_root_entry().get_entry('structure.element_list').set_default_unit('atomic_mass','mass')

        if cartesian:
            cs = self.get_root_entry().get_entry('structure.atom_list.coordinate_system')
            if cs is None:
                cs = PrimitiveEntry('coordinate_system')
                self.get_root_entry().get_entry('structure.atom_list').add_entry(cs)
            cs.set_value('cartesian')
        if angstrom:
            if cartesian:
                self.get_root_entry().get_entry('structure.atom_list.atoms').set_default_unit('angstrom','length')
            self.get_root_entry().get_entry('structure.unit_cell').set_default_unit('angstrom','length')
        else:
            if cartesian:
                self.get_root_entry().get_entry('structure.atom_list.atoms').set_default_unit('bohr','length')
            self.get_root_entry().get_entry('structure.unit_cell').set_default_unit('bohr','length')

        atomtable.reset_data()
        atmblock = self.get_root_entry().get_entry('structure.atom_list.atoms')

        elem_names=[]
        for iatm in range(natms):
            pos = coord[iatm]
            element = symbols[iatm]
            if not element in elem_names:
                elem_names.append(element)
            if not mobile:
                row = {'element':element,'rx':coord_fac * pos[0],'ry':coord_fac * pos[1],'rz':coord_fac * pos[2]}
            else:
                row = {'element':element,'rx':coord_fac * pos[0],'ry':coord_fac * pos[1],'rz':coord_fac * pos[2],
                        'mobile':True}
            atomtable.add_row_data(row)

        a = VectorEntry('a_vector')
        for av in avec:
            a.add_value(av * cell_fac)
        b = VectorEntry('b_vector')
        for bv in bvec:
            b.add_value(bv * cell_fac)
        c = VectorEntry('c_vector')
        for cv in cvec:
            c.add_value(cv * cell_fac)
            blunitcell.add_entry(a)
            blunitcell.add_entry(b)
            blunitcell.add_entry(c)

        elem_table.reset_data()
        for elemname in elem_names:
            if elemname is None or len(elemname)==0:
                continue
            belename = re.sub("\d","",elemname)
            an = atomic_numbers[belename]
            mas = atomic_masses[an]
            if zeta is None:
                row = {'element':elemname,'atomicnumber':an,'mass':mas}
            else:
                if elemname in zeta:
                    row = {'element':elemname,'atomicnumber':an,'mass':mas,'zeta':str(zeta[elemname])}
                else:
                    row = {'element':elemname,'atomicnumber':an,'mass':mas,'zeta':'0.0'}
            elem_table.add_row_data(row)

    def inpfile_exists(self):
        return self.inp_exists

    def get_root_entry(self):
        ''' returns the 'root entry' of the F_INP file. '''
        return self.root_block

    def get_filenames_data(self):
        ''' returns an instance of a file_names.data object '''
        return self.filename_data

    def get_file_name(self):
        ''' returns the file name of F_INP '''
        return self.file_name
    
    def __parse_input(self,inpfile=None):
        inpbuf = self.__get_preprocessed_input(inpfile)
        if inpbuf is None:
            return
        logger.debug('')
        logger.debug("parsing input file : " +self.file_name)
        self.__parse_core(inpbuf)
        logger.debug('parsed input')

    def __parse_core(self,inpbuf):
        parsing_table = False
        table = None
        for li in inpbuf:
            line = li.get_line()
            ent = None
            if line.endswith('{'):
                block = BlockEntry(line.split('{')[0],inp=None,line_no=li.get_line_no())
                self.curr_pblock.add_entry(block)
                self.curr_pblock = block
            elif line=='}':
                self.curr_pblock = self.curr_pblock.get_parent()
                parsing_table = False
                table = None
            elif line.startswith('#units'):
                line = line.replace('#units','')
                ll = line.split()
                for l in ll:
                    self.curr_pblock.set_default_unit(l)
            elif line.startswith('#tag') or line.startswith('#default') or parsing_table:
                parsing_table = True
                if line.startswith('#tag'):
                    if table==None:
                        table = TableEntry(line_no=li.get_line_no())
                    line = line.replace('#tag','').strip()
                    idents = line.split()
                    for ident in idents:
                        table.add_identifier(ident)
                    idents = table.get_identifiers()
                    defs = table.get_defaults()
                    keys = defs.keys()
                    for key in keys:
                        if not key in idents:
                            table.add_identifier(key)
                    self.curr_pblock.add_entry(table)
                elif line.startswith('#default'):
                    if table==None:
                        table = TableEntry(line_no=li.get_line_no())
                    line = line.replace('#default','')
                    defs = line.split(',')
                    for defa in defs:
                        defas = defa.strip().split('=')
                        if len(defas)==2:
                            table.add_default(defas[0].strip(),defas[1].strip())
                else:
                    vals = line.split()
                    iden = table.get_identifiers()
                    row = {}
                    count=0
                    defs = table.get_defaults()
                    for id in iden:
                        if len(vals)>count:
                            if vals[count] != '*':
                                row[id] = vals[count]
                            else:
                                if id in defs:
                                    row[id] = defs[id]
                                else:
                                    row[id] = '*'
                        else:
                            if id in defs:
                                row[id] = defs[id]
                            else:
                                row[id] = '*'
                        count += 1
                    table.add_row_data(row)
            else:
                vals = line.split('=')
                if len(vals)==2:
                    if vals[1].strip().startswith("\"") and vals[1].strip().endswith("\""):
                        ent = PrimitiveEntry(vals[0].strip(),line_no=li.get_line_no())
                        ent.set_value(vals[1].strip()[1:-1])
                        ent.set_double_quote(True)
                    else:
                        va = vals[1].strip().split()
                        if len(va)==1:
                            ent = PrimitiveEntry(vals[0].strip(),line_no=li.get_line_no())
                            ent.set_value(va[0].strip())
                        elif len(va) == 2:
                            ent = PrimitiveEntry(vals[0].strip(),line_no=li.get_line_no())
                            ent.set_value(va[0].strip())
                            ent.set_unit(va[1].strip())
                        elif len(va) >= 3:
                            try:
                                float(va[1].strip())
                                ent = VectorEntry(vals[0].strip(),line_no=li.get_line_no())
                                for v in va:
                                    ent.add_value(v)
                            except ValueError:
                                va[1].strip()
                                ent = PrimitiveEntry(vals[0].strip(),line_no=li.get_line_no())
                                ent.set_value(va[0].strip())
                                ent.set_unit(va[1].strip())
                    if ent != None:
                        self.curr_pblock.add_entry(ent)

    def get_input_as_list(self):
        ''' returns all input entries as a python list. '''
        li = []
        self.__dump(li)
        for l in li:
            logger.debug(l.get_long_name())
        return li

    def __dump(self,li=None):
        if self.root_block==None:
            return None
        entries = self.root_block.get_entries()
        inpstr=''
        for ent in entries:
            if ent.is_block():
                if li!=None:
                    li.append(ent)
                tmpstr=self.__dump_block(ent,li)
                if tmpstr!=None:
                    inpstr += tmpstr
        if li==None:
            logger.debug('input :')
            logger.debug('\n'+inpstr)
        return inpstr

    def __str__(self):
        return self.__dump()

    def __dump_block(self,bloc,li=None):
        inpstr=''
        if li==None:
            logger.debug('long name: '+bloc.get_long_name())
        ents = bloc.get_entries()
        if len(ents)==0:
            return
        depth = bloc.get_block_depth()
        ws = get_ws(depth)
        inpstr += ws+bloc.get_name()+'{\n'
        units = bloc.get_default_units()
        if len(units)!=0:
            inpstr += ws+'  #units'
            for k,v in units.items():
                inpstr += ' '+v
            inpstr += '\n'
        for entry in ents:
            if li!=None:
                li.append(entry)
            if entry.is_block():
                tmpstr=self.__dump_block(entry,li)
                if tmpstr!=None:
                    inpstr += tmpstr
            else:
                if li==None:
                    logger.debug('long name: '+entry.get_long_name())
                entry.set_ws(ws)
                inpstr += str(entry)+'\n'
        inpstr+= ws+'}\n'
        return inpstr

    def __get_preprocessed_input(self,inp=None):
        inpfile = None
        if inp is None and self.file_name is None and not \
            os.path.exists(self.filename_data.get_file_name('F_INP')):
            logger.error('F_INP does not exist.')
        try:
            if self.file_name is None and self.filename_data is not None:
                inpfile = open(self.filename_data.get_file_name('F_INP'))
            elif inp is None:
                inpfile = open(self.file_name)
            elif inp is not None:
                inpfile = open(inp)
        except IOError:
            #logger.error(' input open error.')
            return None
        inpbuf=[]
#        for line in inpfile:
        icount=0
        for ll in inpfile:
            icount = icount+1
            line = Line(ll,icount)
            line.set_line(line.get_line().strip())
            line.set_line(line.get_line().replace('!#','#'))
            if line.get_line().startswith('!') or line.get_line().startswith('//'):
                continue
            lines = line.get_line().split('!')
            lines = lines[0].split('//')
            if(not lines[0].startswith("#")):
                lines = lines[0].split(',')
            for l in lines:
                l=l.strip()
                if len(l)==0:
                    continue
                inpbuf.append(Line(l,line.get_line_no()))
        inpbuf2=[]
        for line in inpbuf:
            line.set_line(line.get_line().strip())
            lines = line.get_line().split('{')
            if len(lines)==1:
                inpbuf2.append(Line(lines[0],line.get_line_no()))
            else:
                endswithlb = line.get_line().endswith('{')
                for ll in lines:
                    ll=ll.strip()
                    if len(ll)==0:
                        continue
                    inpbuf2.append(Line(ll+'{',line.get_line_no()))
                if not endswithlb:
                    inpbuf2[len(inpbuf2)-1].set_line(inpbuf2[len(inpbuf2)-1].get_line()[0:-1])

        inpbuf = []
        for line in inpbuf2:
            line.set_line(line.get_line().strip())
            lines = line.get_line().split('}')
            if len(lines)==1:
                inpbuf.append(Line(lines[0],line.get_line_no()))
            else:
                for i in range(len(lines)-1):
                    ll = lines[i]
                    ll=ll.strip()
                    if len(ll)==0:
                        inpbuf.append(Line('}',line.get_line_no()))
                    else:
                        inpbuf.append(Line(ll,line.get_line_no()))
                        inpbuf.append(Line('}',line.get_line_no()))

        logger.debug('preprocessed input : ')
        for line in inpbuf:
            logger.debug(line)

        inpfile.close()
        return inpbuf

    def save(self,logoffset='',fnames_data=False, ppdict=None):
        ''' save the input to filename.
            if direc is given, will save to the corresponding file under
            that directory'''

        fname = self.get_filenames_data().get_file_name('F_INP')
        if fname is None:
            fname = self.get_file_name()
        logger.debug(logoffset+'dumping F_INP to : '+fname)
        try:
            fi = open(fname,'w')
            fi.write(str(self))
            fi.flush()
            fi.close()
        except IOError:
            logger.error(logoffset+'failed write to: '+fname)

        if fnames_data:
            elems = self.get_defined_elements()
            for i,elem in enumerate(elems):
                if ppdict is not None and elem in ppdict:
                    self.filename_data.set_file_name('F_POT('+str(i+1)+')',ppdict[elem])
                else:
                    belename = re.sub("\d","",elem)
                    self.filename_data.set_file_name('F_POT('+str(i+1)+')', 
                                                 get_default_ppfile_from_dir(belename,full_path=True))
            self.filename_data.save()

    def get_defined_elements(self):
        elem_table = self.get_root_entry().get_entry('structure.element_list.table')
        if elem_table is None:
            logger.error('element table undefined')
            return None
        nelem = elem_table.get_num_data()
        ret = []
        for i in range(nelem):
            ename = elem_table.get_data(i,'element')
            ret.append(ename)
        return ret

class InputEntry:
    ''' the base class for input entries '''
    def __init__(self,name,line_no=-1):
        ''' parameters : name the name for this entry '''
        self.name = name
        self.type = -1
        self.ws = ''
        self.parent=None
        self.spec=None
        self.line_no=line_no

    def set_line_no(self,line_no):
        ''' set the line no associated to this input entry.'''
        self.line_no=line_no

    def get_line_no(self):
        ''' get the line no associated to this input entry.'''
        return self.line_no

    def get_name(self):
        ''' returns the name of this entry '''
        return self.name

    def set_specification(self,spe):
        ''' set the dictionary object which specifies
            the present input entry. '''
        self.spec = spe

    def get_specification(self):
        ''' returns the specification dictionary.
            will return None if no specification is set. '''
        return self.spec

    def get_type(self):
        ''' returns the 'type' of this entry;
            it will be one of BLOCK, TABLE, PRIMITIVE or VECTOR
            subclass MUST implement this method. '''
        raise NotImplementedError()

    def get_parent(self):
        ''' get the reference to the parent block '''
        return self.parent

    def _set_parent(self,parent):
        self.parent = parent

    def get_long_name(self):
        ''' get the 'long name' of this input entry.
            examples of a 'long name' :
            control.cpumax
            accuracy.ksampling
            accuracy.ksampling.method
            structure.atom_list.atoms.table
            ...
            ... '''
        ret=[]
        ent=self
        while 1:
            if ent==None:
                break
            ret.append(ent.get_name())
            ent = ent.get_parent()
        ret.reverse()
        retret=''
        for re in ret:
            retret += re+'.'
        retret = retret[0:-1].replace('root.','').lower()
        return retret

    def is_block(self):
        ''' returns True if this input entry is of type BLOCK '''
        return self.get_type()==BLOCK

    def set_ws(self,ws):
        ''' sets the amount of white space for printing
            parameter : ws the white-space string '''
        self.ws = ws

    def delete_me(self):
        ''' remove this entry from its parent block
            (this method will not work for the root block,
            whose parent is None). '''
        if self.parent != None:
            self.parent.remove_entry(self)

class BlockEntry(InputEntry):
    ''' a class which represents a BLOCK entry '''
    def __init__(self,name,inp=None,line_no=-1):
        InputEntry.__init__(self,name,line_no)
        self.default_unit = {}
        self.entries = []
        self.input=None
        if inp != None:
            self.input = inp

    def get_input_obj(self):
        if self.input!=None:
            return self.input
        ent=self
        while 1:
            if ent==None:
                return None
            if ent.input!=None:
                return ent.input
            ent = ent.get_parent()

    def get_type(self):
        return BLOCK

    def set_default_units(self,du):
        ''' set the 'default unit' associated with this block.
            parameter: du a Python dictionary;
            example: du[unit_type_force] >> 'hartree/bohr' '''
        self.default_unit = du

    def set_default_unit(self,uval,typ=None):
        ''' set the 'default_unit' the unit type will be
            automatically resolved unless you give it
            explicitly.
            parameter : uval the value of the unit
            typ : the unit type (optional) '''
        if typ==None:
            typ = resolve_unit(uval)
            if typ != None:
                self.default_unit[typ]=str(uval)
        else:
            self.default_unit[typ] = str(uval)

    def get_default_unit(self,utyp):
        ''' get the 'default unit'
            will return None if the specified unit-type
            does not exist.
            parameter : utyp the type of the unit to be obtained. '''
        if not utyp in self.default_unit:
            return None
        return self.default_unit[utyp]

    def get_default_units(self):
        ''' get the 'default unit' dictionary itself. '''
        return self.default_unit

    def add_entry(self,entry):
        ''' add a new entry to this block. if an entry of the same name
            already exist, it will be replaced.
            parameter : entry the entry to be added.
            it must be an instance of a derived class of InputEntry. '''
        if not isinstance(entry,InputEntry):
            logger.error(type(entry))
            raise TypeError('only an instance of InputEntry can be added.')

        ents = self.get_entries()
        nam=entry.get_name().lower()
        exi=False
        for ent in ents:
            if ent.get_name().lower()==nam:
                self.entries[self.entries.index(ent)]=entry
                exi=True
                break
        if not exi:
            self.entries.append(entry)
        entry._set_parent(self)
        lname = entry.get_long_name()
        inp=self.get_input_obj()

    def remove_entry(self,entry):
        ''' remove the specified entry (if it exists)
            parameter : entry the entry to be removed. '''
        if entry in self.entries:
            self.entries.remove(entry)

    def get_entries(self):
        ''' get all the entries contained in this block. '''
        return self.entries

    def get_entry(self,entry_name):
        ''' get the entry
            parameter : entry_name the name of the entry to be obtained;
            (case insensitive). if entry_name is delimited by a '.', then
            get_entry will recursively call itself assuming that the list obtained by
            entry_name.split('.')[:1] represent a block of the corresponding depth.
            will return None if no such entry exists. '''
        entnam = entry_name.split(".")
        if len(entnam)<=1:
            for ent in self.entries:
                if ent.get_name().lower()==entry_name.lower():
                    return ent
        else:
            ent = self.get_entry(entnam[0])
            if ent == None:
                return None
            str = ".".join(entnam[1:])
            return ent.get_entry(str)
        return None

    def entry_count(self):
        ''' get the number of entries defined in this block '''
        return len(self.entries)

    def get_block_depth(self):
        ''' get the 'depth' of the block '''
        parent = self
        depth = 0
        while 1:
            parent = parent.get_parent()
            if parent == None:
                return depth
            depth += 1
    def __str__(self):
        return self.get_name()

class PrimitiveEntry(InputEntry):
    ''' this class represents a 'primitive entry'.
        a primitive entry is specified by
        name = value unit
        a value can be an int, float, bool or string;
        only a float can have units. '''

    def __init__(self,name,val='',uni='',dq=False,line_no=-1):
        ''' parameter : name the name of this entry '''
        InputEntry.__init__(self,name,line_no)
        self.unit = str(uni)
        self.value = str(val)
        self.dq = dq

    def get_type(self):
        return PRIMITIVE

    def set_unit(self,unit):
        ''' set the unit for this entry
            paremeter : unit the unit; it must be a string '''
        self.unit = str(unit)

    def set_value(self,value):
        ''' set the value for this entry
            parameter : value the value; it must be a string. '''
        tmpstr=str(value)
        if str(value)=='False':
            tmpstr='off'
        elif str(value)=='True':
            tmpstr='on'
        self.value = tmpstr

    def get_unit(self):
        ''' get the unit '''
        return self.unit

    def get_value(self):
        ''' get the value. if the specification is given,
            will return a suitably converted value. '''
        return get_value(self.spec,self.value)

    def get_str_value(self):
        ''' get the value in string. '''
        return str(self.value)

    def set_double_quote(self,dq):
        ''' if set to True, will double quote the value. '''
        self.dq = dq

    def is_double_quoted(self):
        ''' returns True if the value should be double-quoted on save. '''
        return self.dq

    def __str__(self):
        if not self.dq:
            return self.ws+'  '+self.name+' = '+str(self.value)+' '+self.unit
        else:
            return self.ws+'  '+self.name+' = \"'+str(self.value)+'\"'

class VectorEntry(InputEntry):
    ''' this class represents a 'vector entry'. note that
        the only usage of this data structure is for the unit cell. '''
    def __init__(self,name,line_no=-1):
        InputEntry.__init__(self,name,line_no)
        self.value=[]

    def set_value(self,value):
        ''' set the value.
            parameter : value the value of
            the vector entry. it must be a Python list. '''
        self.value = value

    def add_value(self,va):
        ''' add values to the vector
            parameter : va the value to be added '''
        self.value.append(str(va))

    def get_value(self):
        ''' get the value.
            returns a python list. '''
        if self.spec==None:
            return self.value
        ret=[]
        inde=0
        for val in self.value:
            ret.append(get_value(self.spec,val,inde))
            inde=inde+1
        return ret

    def get_type(self):
        return VECTOR

    def __str__(self):
        ret = self.name+' ='
        for val in self.value:
            ret += ' '+str(val)
        return self.ws+'  '+ret

class TableEntry(InputEntry):
    ''' this class represents a tabular data.
        refer to the PHASE manual for details of
        this data atructure'''
    def __init__(self,line_no=-1):
        InputEntry.__init__(self,'table',line_no)
        self.identifiers = []
        self.data = []
        self.default = {}

    def reset_data(self):
        self.data = []

    def get_num_data(self):
        ''' returns the number of row data included in this table '''
        return len(self.data)

    def get_type(self):
        return TABLE

    def add_default(self,iden,val):
        ''' a table can have 'default values'.
            the default value will be stored in a Python dictionary
            parameter : iden the identifier for the default value.
            parameter : val  the actual value of the default value. '''
        logger.debug(iden)
        logger.debug(val)
        self.default[iden] = str(val)

    def get_defaults(self):
        return self.default

    def add_identifier(self,ident):
        ''' add table identifiers
            (note that the table identifiers are stored in lists). '''
        self.identifiers.append(ident)

    def set_identifiers(self,ident):
        self.identifiers = ident

    def get_identifiers(self):
        ''' returns a list of column identifiers. '''
        return self.identifiers

    def valid_identifier(self,iden):
        ''' return True if iden is included in the list of identifiers, False otherwise '''
        return iden in self.identifiers

    def get_table_data(self):
        ''' returns the table data it self.
            table data are stored in a list of Python dictionaries.
            for instance, get_table_data()[0]['element'] should return
            the element of the 0-th atom. '''
        if self.spec==None:
            return self.data
        ret=[]
        for k in range(len(self.data)):
            self.data[k] = self.get_row_data(k)
        return self.data

    def get_str_table_data(self):
        return self.data

    def get_corresponding_data(self,key,val,ident):
        ''' similar to get_data, except that the row index is resolved differently.
            the row whose 'key' column equals 'val' will be adopted.
            will return None if 'key' is unfound in the identifiers,
            or if no  equals 'val' '''
        if not key in self.identifiers:
            return None
        for i in range(len(self.data)):
            rowdata = self.data[i]
            if not key in rowdata.keys():
                return None
            coldat = rowdata[key]
            if coldat == val:
                return get_data(i,ident)
        return None
        
    def get_data(self,row,ident):
        ''' returns the specified data
            parameter : row the row index of the data.
            ident : the identifier for the data.
            if the specification file is given, a suitably converted
            value will be returned. '''

        if row>=len(self.data):
            return None
        if not ident in self.data[row]:
            return None

        if self.spec==None:
            return self.data[row][ident]
        cols=self.spec[columns].split(',')
        if not ident.lower() in cols:
            return self.data[row][ident]
        inde = cols.index(ident)
        return get_value(self.spec,self.data[row][ident],inde)

    def get_row_data(self,row):
        ''' get the specified row-data
            parameter : row the row-index of the data
            if the specification file is given, a suitably converted
            value will be returned.'''
        if row>=len(self.data):
            return None

        if self.spec==None:
            return self.data[row]

        keys=self.data[row].keys()
        for key in keys:
            self.data[row][key]=self.get_data(row,key)
        return self.data[row]

    def get_str_row_data(self,row):
        ''' get the row data as it is.'''
        if row>=len(self.data):
            return None
        return self.data[row]

    def add_row_data(self,row):
        ''' add row data.
            remark: dat must be a python dictionary '''
        rkeys = row.keys()
        for rkey in rkeys:
            row[rkey] = str(row[rkey])

        for rkey in rkeys:
            if not rkey in self.identifiers:
                self.identifiers.append(rkey)
        self.data.append(row)

    def remove_row_data(self,rowdata):
        ''' remove the specified row. will do nothing if the
            specified row-data does not exist.
            parameter : the row data to be removed. '''
        if rowdata in self.data:
            self.data.remove(rowdata)

    def remove_row_data_by_index(self,index):
        if index>=len(self.data):
            return
        self.data.pop(index)

    def set_data(self,row,ident,value):
        if row>=len(self.data):
            logger.warn('at set_data : row index is larger than data.')
            logger.warn('nothing to do')
            return
        self.data[row][ident] = str(value)

    def set_row_data(self,row,dat):
        ''' remark: dat must be a python dictionary '''
        if row>=len(self.data):
            logger.warn('at set_data : row index is larger than data.')
            logger.warn('nothing to do')
            return
        for da in dat:
            dat[da] = str(dat[da])
        self.data[row] = dat

    def __str__(self):
        ret = self.ws+'  #tag'
        for iden in self.identifiers:
            ret += '    '+str(iden).center(17)
        for row in self.data:
            ret += '\n'+self.ws+'    '
            for iden in self.identifiers:
                if iden in row:
                    tmpstr=str(row[iden])
                    if str(row[iden])=='False':
                        tmpstr='off'
                    elif str(row[iden])=='True':
                        tmpstr='on'
                    ret += '    '+tmpstr.rjust(17)
                else:
                    ret += '    *'
        return ret

class FileNamesData:
    fnames_keys = [
     'F_INP'
    ,'F_POT(1)'
    ,'F_POT(2)'
    ,'F_POT(3)'
    ,'F_POT(4)'
    ,'F_POT(5)'
    ,'F_POT(6)'
    ,'F_POT(7)'
    ,'F_POT(8)'
    ,'F_POT(9)'
    ,'F_POT(10)'
    ,'F_POT(11)'
    ,'F_POT(12)'
    ,'F_POT(13)'
    ,'F_POT(14)'
    ,'F_POT(15)'
    ,'F_POT(16)'
    ,'F_CHR'
    ,'F_WFk'
    ,'F_WANNIER'
    ]

    neb_file_idents={
      'F_NEB_OUT':'output_neb'
     ,'F_IMAGE(0)':'endpoint0.data'
     ,'F_IMAGE(-1)':'endpoint1.data'
     ,'F_NEB_STOP':'./nfnebstop.data'
     ,'F_NEB_CNTN':'./neb_continue.data'
     ,'F_NEB_ENF':'./nfnebenf.data'
     ,'F_NEB_DYNM':'./nfnebdynm.data'
    }

    fname_idents={
     'F_INP':'./nfinp.data'
    ,'F_POT(1)':''
    ,'F_POT(2)':''
    ,'F_POT(3)':''
    ,'F_POT(4)':''
    ,'F_POT(5)':''
    ,'F_POT(6)':''
    ,'F_POT(7)':''
    ,'F_POT(8)':''
    ,'F_POT(9)':''
    ,'F_POT(10)':''
    ,'F_POT(11)':''
    ,'F_POT(12)':''
    ,'F_POT(13)':''
    ,'F_POT(14)':''
    ,'F_POT(15)':''
    ,'F_POT(16)':''
    ,'F_PKB':'./vkb.data'
    ,'F_PD':'./vd.data'
    ,'F_PPC':'./vpc.data'
    ,'F_STOP':'./nfstop.data'
    ,'F_CNST':'./nfcnst.data'
    ,'F_OPGR':'./opgr.data'
    ,'F_MATBP':'./matrix.BP'
    ,'F_KPOINT':'./kpoint.data'
    ,'F_KINDEX':'./f.kp0'
    ,'F_SPG':'./bnprp4.i5'
    ,'F_KPGN':'./bnkpgn.i5'
    ,'F_OTP':'./nfotp.data'
    ,'F_DYNM':'./nfdynm.data'
    ,'F_EGRID':'./nfegrid.data'
    ,'F_LDOS':'./nfldos.data'
    ,'F_ENF':'./nfefn.data'
    ,'F_GPT':'./nfgpt.data'
    ,'F_CHGT':'./nfchgt.data'
    ,'F_CHGO':'./nfchgo.data'
    ,'F_CHGU':'./nfchgu.data'
    ,'F_ENERG':'./nfenergy.data'
    ,'F_CNTN':'./continue.data'
    ,'F_CNTN_BIN':'./continue_bin.data'
    ,'F_ZAJ':'./zaj.data'
    ,'F_VLC':'./nfvlc.data'
    ,'F_CNTN_BIN_STM':'./continue_bin_stm.data'
    ,'F_DOS':'./dos.data'
    ,'F_CHR':'./nfchr.cube'
    ,'F_WFk':'./nfwfk.cube'
    ,'F_WANNIER':'./nfwannier.cube'
    ,'F_CNTN_WAN':'./nfcontinue_wannier.data'
    ,'F_POT_WAN':'./nfpotential_wannier.data'
    ,'F_ELF':'./nfelf.data'
    ,'F_BERRY':'./berry.data'
    ,'F_EFFCHG':'./effchg.data'
    ,'F_FORCE':'./force.data'
    ,'F_MODE':'./mode.data'
    ,'F_EPSILON':'./epsilon.data'
    ,'F_EPSOUT':'./eps.data'
    ,'F_NLO':'./nlo.data'
    ,'F_MAGOPT':'./magopt.data'
    ,'F_PHDOS':'./phdos.data'
    ,'F_PHDOS':'./phdos.data'
    ,'F_OCCMAT':'./occmat.data'
    ,'F_STRFRC':'./strfrc.data'
    ,'F_PSTRN':'./positron.cube'
    ,'F_VELEC':'./electron.cube'
    ,'F_EPPAIR':'./ep_pair.cube'
    ,'F_VELEC_GRAD':'./electron_grad.cube'
    ,'F_EFERMI':'./nfefermi.data'

    ,'F_ZAJ_filetype':'unified'
    ,'F_CHGT_filetype':'unified'
    ,'F_CNTN_BIN_filetype':'unified'
    ,'F_CNTN_filetype':'unified'

    ,'F_ZAJ_in_filetype':'nogiven'
    ,'F_CHGT_in_filetype':'nogiven'
    ,'F_CNTN_BIN_in_filetype':'nogiven'
    ,'F_CNTN_in_filetype':'nogiven'

    ,'F_CNTN_BIN_PAW':'./continue_bin_paw.data'
    ,'F_CNTN_BIN_RSNLPP':'./continue_bin_rsnlpp.data'
    }

    mpifiletypes_idents = {
    'F_ZAJ_filetype' : 'unified'
    ,'F_CHGT_filetype' : 'unified'
    ,'F_CNTN_BIN_filetype' : 'unified'
    ,'F_CNTN_filetype' : 'unified'

    ,'F_ZAJ_in'      : './zaj.data'
    ,'F_CNTN_in'     : './continue.data'
    ,'F_CNTN_BIN_in' : './continue_bin.data'
    ,'F_CHGT_in'     : './nfchgt.data'

    ,'F_ZAJ_bak'      : './zaj.data'
    ,'F_CNTN_bak'     : './continue.data'
    ,'F_CNTN_BIN_bak' : './continue_bin.data'
    ,'F_CHGT_bak'     : './nfchgt.data'

    ,'F_ZAJ_in_filetype' : 'nogiven'
    ,'F_CHGT_in_filetype' : 'nogiven'
    ,'F_CNTN_BIN_in_filetype' : 'nogiven'
    ,'F_CNTN_in_filetype' : 'nogiven'
    }

    f_param_name_idents = {
    'F_PRM':'./m_ArraySize_Parameters.f90'
    }

    default_identifiers ={
      '&fnames':fname_idents
     ,'&nebfiles':neb_file_idents
     ,'&mpifiletypes':mpifiletypes_idents
     ,'&f_param_name':f_param_name_idents
    }

    def __init__(self,inp=None,nfinp=None):
        self.fnames_data_exists=False
        self.inp = inp
        self.file_name = join('.','file_names.data')
        if not os.path.exists(self.file_name):
            self.fnames_data_exists=False
        else:
            self.fnames_data_exists=True
            #return
            #if nfinp is None:
            #    logger.debug('file_names.data does not exist.')
            #    self.valid_fnames=False
            #    return
            #else:
            #    f=open(self.file_name,'w')
            #    f.write('&fnames\n')
            #    f.write("F_INP='"+nfinp+"'\n")
            #    f.write("/")
            #    f.flush()
            #    f.close()
        self.entry = {}
        self.entry["&fnames"]=self.default_identifiers
        if self.fnames_data_exists:
            self.valid_fnames = self._ini()

    def is_valid(self):
        return self.valid_fnames

    def get_file_name(self,ident,file_set="&fnames",to_upper=True):
        ''' get the file name corresponding to the specified identifier.
            parameter : ident the identifier '''
        ide = str(ident)
        try:
            if to_upper:
                ide = str(ident).upper()
            return self.entry[file_set][ide]
        except:
            return None

    def set_file_name(self,ident,file_name,file_set="&fnames"):
        ''' set the file name corresponding to the specified identifier.
            parameter : ident the identifier.
            parameter : file_name the filename for the file corresponding to ident. '''
        f = self.entry[file_set]
        f[ident] = file_name
        try:
            logger.debug(self.fnames_keys.index(ident))
        except ValueError:
            logger.debug("added "+ident+" to the list of keys")
            self.fnames_keys.append(ident)

    def save(self,logoffset=''):
        ''' save the file_names.data file to disk.
            if direc is given, the method will save to the file_names.data
            under dir. '''
#        fname = join(self.inp.get_curr_dir(),'file_names.data')
        fname = 'file_names.data'

        logger.debug(logoffset+"dumping file_names.data to : "+fname)

        try:
            f=open(fname,"w")
            string = self.__str__()
            f.write(string)
            f.flush()
            f.close()
            ekeys = self.entry.keys()
        except IOError:
            logger.error(logoffset+"failed write to "+self.file_name)

    def __str__(self):
        ret=''
        ekeys = self.entry.keys()
        for ekey in ekeys:
            ret += ekey+'\n'
            iden = self.entry[ekey]
            if ekey=='&fnames':
                for key in self.fnames_keys:
                    if key in iden and len(iden[key])!=0:
                        ret += key+' = \''+iden[key]+'\'\n'
            else:
                ikeys = iden.keys()
                for ikey in ikeys:
                    if ikey in iden:
                        ret += ikey+' = \''+iden[ikey]+'\'\n'
            ret += '/\n'
        return ret

    def file_names_data_exists(self):
        return self.fnames_data_exists

    def _ini(self):
        logger.debug("")
        logger.debug("checking the file_names.data file...")
        boo = True
        fi = open(self.file_name,"r")
        currident = ''
        for line in fi:
            line = line.strip()
            if line.startswith("&"):
                if line in self.default_identifiers.keys():
                    self.entry[line] = self.default_identifiers[line]
                else:
                    self.entry[line] = {}
                currident = line
            elif line.startswith("/"):
                currident = ''
            else:
                if not len(currident)==0:
                    bb = line.split("=")
                    ident = ''
                    val = ''
                    if len(bb)>1:
                        ident = bb[0].strip()
                        if len(bb[1])<1:
                            continue
                        if ident in self.entry[currident].keys():
                            val = bb[1].strip()[1:len(val)-1]
                            self.entry[currident][ident] = val
                        else:
                            logger.error("  unknown file pointer : "+ident.rjust(15)+" in "+currident)
                            boo=False
        if boo:
            logger.debug("no errors were found in the file_names.data file")
        else:
            logger.error("encountered errors in the file_names.data file")

        return boo

    def debug_output(self):
        logger.debug(self.file_name)
        ekeys = self.entry.keys()
        for ekey in ekeys:
            logger.debug('entry : '+ekey)
            iden = self.entry[ekey]
            if ekey=='&fnames':
                for key in self.fnames_keys:
                    if key in iden:
                        logger.debug('  '+key.ljust(10)+' : '+iden[key])
            else:
                ikeys = iden.keys()
                for ikey in ikeys:
                    if ikey in iden:
                        logger.debug('  '+ikey.ljust(14)+' : '+iden[ikey])

    def convert_pp_to_abspath(self):
        ''' convert the path to the pseudo-potential file to absolute'''
        if not '&fnames' in self.entry:
            logger.warn(' no &fname entry in file_names.data')
            return
        fnams = self.entry['&fnames']
        idents = fnams.keys()
        for ident in idents:
            if ident.startswith('F_POT('):
                fnams[ident] = os.path.abspath(fnams[ident])

    def convert_to_abspath(self,ident):
        ''' convert the path to file specified by ident. '''
        if not '&fnames' in self.entry:
            logger.warn(' no &fname entry in file_names.data')
            return
        fnams = self.entry['&fnames']
        if not ident in fnams:
            logger.warn('could not find '+ident+' in file_names.data')
            return
        fnams[ident] = os.path.abspath(fnams[ident])
