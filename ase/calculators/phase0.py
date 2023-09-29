"""This module defines an ASE interface to PHASE/0.

https://azuma.nims.go.jp/cms/
"""

import sys,os
from ase.io.phase0 import Input,BlockEntry,PrimitiveEntry,get_pp_path_as_list, \
    read_phase0_out, update_continue_data, get_outputxxx, get_jobstatusxxx, restartable, \
    get_mostrescent_outputxxx,Hartree_Bohr3_2_eV_Ang3
from ase.calculators.calculator import FileIOCalculator, Parameters, kpts2mp, \
    ReadError
import numpy as np


class phase0(FileIOCalculator):
    """Class for doing PHASE/0 calculations."""
    implemented_properties = ['energy', 'forces', 'stress']
    #command = 'phase'
    default_parameters = dict(
        ecut='25',
        ecut_charge=None,
        xc='gga',
        smearing=None,
        symmetry=True,
        inversion=False,
        stress=False,
        stress_correction=True,
        kpts=3.5,
        gamma=False,
        mobile=False,
        spin=False,
        zeta=None,
        preproc=False,
        dE=1.0e-9,
        pp=None)
    phaseinp = None

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='phase0', atoms=None, scratch=None, from_scratch=True, **kwargs):
        """Construct PHASE/0-calculator object.  """

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        self.__init__phase(from_scratch)
        self.phaseinp = self.get_input()

    def clean(self):
        remtarget = ['nfefn.data___','nfdynm.data___','zaj.data','nfchgt.data','continue.data',
                     'continue_bin.data','continue_bin_paw.data']
        for remt in remtarget:
            if os.path.exists(remt):
                os.remove(remt)
        jobs = get_jobstatusxxx()
        outs = get_outputxxx()
        for job in jobs:
            os.remove(job)
        for out in outs:
            os.remove(out)
        
    def __init__phase(self,from_scratch=True):
        if from_scratch:
            self.clean()

    def export_input(self, atoms, properties=None, system_changes=None, default_fnames=True):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        if system_changes is None:
            system_changes = []
        finp = self.directory+'/nfinp.data'
        if not default_fnames:
            finp = self.directory+'/nfinp.data___'
        param = self.parameters
        if self.phaseinp is None:
            self.phaseinp = Input(coord_file=finp)

        root = self.phaseinp.get_root_entry()
        control = root.get_entry('control')
        if control is None:
            control = BlockEntry('control')
            root.add_entry(control)

        if restartable() and not 'cell' in system_changes:
            control.add_entry(PrimitiveEntry('condition','automatic'))
            self.phaseinp.save(fnames_data=True,ppdict=param.get('pp'))
            update_continue_data(atoms)
            return

        if restartable() and 'cell' in system_changes:
            control.add_entry(PrimitiveEntry('condition','initial'))
        else:
            control.add_entry(PrimitiveEntry('condition','automatic'))
        if param.get('preproc'):
            control.get_entry('condition').set_value('preparation')
        cpum = param.get('cpumax')
        if cpum is not None:
            control.add_entry(PrimitiveEntry('cpumax',cpum,'hour'))
        self.phaseinp.swap_atomic_configuration(atoms,to_file=finp,cartesian=False,
                                                angstrom=True,mobile=param.get('mobile'),zeta=param.get('zeta'))
        blaccuracy = root.get_entry('accuracy')
        if blaccuracy is None:
            blaccuracy = BlockEntry('accuracy')
            root.add_entry(blaccuracy)
        blaccuracy.add_entry(PrimitiveEntry('cutoff_wf',val=str(param.get('ecut')),uni='rydberg'))
        fac = 4.0
        if self.contains_uspp(self.phaseinp):
            fac = 9.0
        blaccuracy.add_entry(PrimitiveEntry('cutoff_cd',val=str(float(param.get('ecut'))*fac),uni='rydberg'))
        if param.get('ecut_charge') is not None:
            blaccuracy.get_entry('cutoff_cd').set_value(param.get('ecut_charge','rydberg'))

        conver = blaccuracy.get_entry('scf_convergence')
        if conver is None:
            conver = BlockEntry('scf_convergence')
            blaccuracy.add_entry(conver)
        conver.add_entry(PrimitiveEntry('delta_total_energy',str(param.get('dE'))))

        structure = root.get_entry('structure')
        symm = structure.get_entry('symmetry')
        if symm is None:
            symm = BlockEntry('symmetry')
            structure.add_entry(symm)
        if param.get('symmetry'):
            symm.add_entry(PrimitiveEntry('method','automatic'))
        else:
            symm.add_entry(PrimitiveEntry('method','manual'))
        if param.get('inversion'):
            symm.add_entry(PrimitiveEntry('sw_inversion','on'))
        if param.get('symmetry_check_criterion'):
            symm.add_entry(PrimitiveEntry('symmetry_check_criterion',str(param.get('symmetry_check_criterion'))))

        if param.get('spin'):
            structure = root.get_entry('structure')
            structure.add_entry(PrimitiveEntry('magnetic_state','ferro'))

        strevl = root.get_entry('structure_evolution')
        if strevl is None:
            strevl = BlockEntry('structure_evolution')
            root.add_entry(strevl)
        stressb = BlockEntry('stress')
        strevl.add_entry(stressb)
        stress = PrimitiveEntry('sw_stress')
        stressb.add_entry(stress)
        if param.get('stress'):
            stress.set_value('on')
            if param.get('stress_correction'):
                stress_c = PrimitiveEntry('sw_smear_ke','on')
                stressb.add_entry(stress_c)
        else:
            stress.set_value('off')

        ksamp = BlockEntry('ksampling')
        blaccuracy.add_entry(ksamp)

        if param.get('gamma'):
            ksamp.add_entry(PrimitiveEntry('method','gamma'))
        else:
            mesh = BlockEntry('mesh')
            ksamp.add_entry(mesh)
       
            mp = None
            if param.kpts is None or isinstance(param.kpts,(float,int)):
                mp = kpts2mp(atoms, param.kpts)
            else:
                mp = param.kpts
            nx = PrimitiveEntry('nx',str(mp[0]))
            ny = PrimitiveEntry('ny',str(mp[1]))
            nz = PrimitiveEntry('nz',str(mp[2]))
            mesh.add_entry(nx)
            mesh.add_entry(ny)
            mesh.add_entry(nz)

        initial_chg = PrimitiveEntry('initial_charge_density','atomic_charge_density')
        blaccuracy.add_entry(initial_chg)
        if not default_fnames:
            fdata = self.phaseinp.get_filenames_data()
            fdata.set_file_name('F_INP','nfinp.data___')
            fdata.set_file_name('F_ENF','nfefn.data___')
            fdata.set_file_name('F_DYNM','nfdynm.data___')
        self.phaseinp.save(fnames_data=True,ppdict=param.get('pp'))

    def write_input(self, atoms, properties=None, system_changes=None, default_fnames=False):
        self.export_input(atoms,properties,system_changes,default_fnames=default_fnames)

    def read_results(self):
        atoms = read_phase0_out(index=-1,fname_from_fnamesdata=True)
        self.results['energy'] = atoms.get_potential_energy()
        forces = atoms.get_forces()
        self.results['forces'] = np.array(atoms.get_forces())
        if self.parameters.get('stress'):
            output = get_mostrescent_outputxxx()
            f = open(output)
            in_stress = False
            strcount = 0
            stress = []
            for line in f:
                line = line.strip()
                if line.startswith('Total STRESS TENSOR'):
                    in_stress = True
                    continue
                if in_stress:
                    words = line.split()
                    stress.append(np.array([-float(words[0])*Hartree_Bohr3_2_eV_Ang3,
                                            -float(words[1])*Hartree_Bohr3_2_eV_Ang3,
                                            -float(words[2])*Hartree_Bohr3_2_eV_Ang3]))
                    if len(stress) == 3:
                        in_stress = False
                        continue
            f.close()
            self.results['stress'] = np.array(stress)

    def contains_uspp(self,input):
        pps = get_pp_path_as_list(input)
        for pp in pps:
            if len(pp.split('_us_'))>1:
                return True
        return False

    def get_input(self):
        if self.phaseinp is None:
            finp = self.directory+'/nfinp.data'
            self.phaseinp = Input(coord_file=finp)
        return self.phaseinp

