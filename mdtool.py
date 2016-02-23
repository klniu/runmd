#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
'''
@Purpose:  Generate command file, topology file(.top), mdp file(.mdp) which are needed by running molecular simulation of surfactant by gromacs.
        -t, --top       Generate topology file(.top)
        -c, --command   Generate command file(.top)
        -m, --mdp       Generate mdp file(.mdp)
'''

import configparser
import argparse
import os.path
from collections import OrderedDict


class ConfigReader(configparser.ConfigParser):
    '''
    Read configure.ini file and initialize the parameters of md.

    The ini file includes several parts:

    [general]: The parameters related the machine and runtime.
    [input]: The start parameters, including forcefield, box and some file.
    [component]: About the component substance in md system.
    [ions]: About the ions.
    [pr]: About the position restraints simulation.
    [md]: The md parameters, including the sequence of md course.
    [em_default]: The default mdp values for energy minimuze, used when no assigned for these values.
    [md_default]: The default mdp values for molecular simulation, used when no assigned for these values.
    [Other Section]: Other Section are em or md sequence defined in [md]. The parameters included below are em or md parameters.

    Please read configure.ini for more information.

    This class will supply these members:
        self.secs: contains sections in ini file except the simulation sections which are is [md] ems and mds keys&values. e.g.:
        self.secs(dict):
            general(dict):
                ...
            input(dict):
                ...
            component(dict):
                ...
            ions(dict):
                ...
            pr(dict):
                ...
            md(dict):
                ...
            analysis(dict):
                ...
        self.ems(OrderedDict): contains simulation sections in [md] ems keys&values.
            emsteep/emcg(dict):
                type:steep/cg
                [nodes]: int
                [mdp]:./include/[section].mdp
                [analysis_dir]:/result/[section]
                ...
        self.mds(OrderedDict): contains simulation sections in [md] mds keys&values.
            nvt(dict):
                type:npt
                [nodes]: int
                [mdp]:./include/[section].mdp
                [analysis_dir]:/result/[section]

                nsteps
                continuation
                tcoupl
                tau-t
                ref-t
                gen_vel
                [gen_temp]
                ...
            npt(dict):
                nvt(dict) +
                pcoupl
                tau-p
                ref-p
                compressibility
                ...

            nvt_pr(dict):
                nvt(dict) +
                pr_resi

                define

            npt_pr(dict):
                npt(dict) +
                pr_resi

                define
    '''

    def __init__(self, f):
        super(ConfigReader, self).__init__(inline_comment_prefixes=(';'))

        try:
            self.read_file(f)
        except configparser.DuplicateSectionError:
            fail('Error: Section Duplicate, please check ini file.')
        except configparser.DuplicateOptionError:
            fail('Error: Key DuplicateKey, please check ini file.' + configparser.DuplicateOptionError.args)
        except:
            fail('Read configuration file Error.')
        self.__check_ini()

        # get keys and values from every secion
        secs = ('general', 'input', 'component', 'ions', 'pr', 'md', 'analysis')
        self.secs = {}
        for sec in secs:
            self.secs[sec] = {k: v for (k, v) in self[sec].items()}

        self.__init_general()
        self.__init_input()
        self.__init_component()
        self.__init_ions()
        self.__init_pr()
        self.__init_md()
        self.__init_analysis()

        # simulations
        self.ems = OrderedDict()
        self.mds = OrderedDict()
        # Initialize energy minumize course
        for md in self.secs['md']['ems'] + self.secs['md']['mds']:
            if not self.has_section(md):
                fail('There is no [{0}] section in ini file. This section is a course in [md] ems values.'.format(md))
            self.__init_simulation(md)

    def __check_ini(self):
        '''Check the ini file'''
        # Default sections and keys in every section.
        secs = {}
        secs['general'] = ('precision', 'nodes', 'opt_mpi', 'maxwarn')
        secs['input'] = ('forcefield', 'start_pdb', 'top', 'box')
        secs['component'] = ('itps', 'mols_num', 'mols_name', 'extra_itps')
        secs['ions'] = ('name', 'num')
        secs['pr'] = ('residures', 'itps')
        secs['md'] = ('ems', 'mds')
        secs['analysis'] = ('rdf',)

        #检查ini配置文件是否有上述节与键
        for sec in secs.keys():
            if not self.has_section(sec):
                fail('There is no {0} section in ini file.'.format(sec))
            for key in secs[sec]:
                if not self.has_option(sec, key):
                    fail('There is no {0} key in {1} section.'.format(key, sec))

    def __init_general(self):
        '''Initialize keys in [general]'''
        try:
            self.secs['general']['nodes'] = int(self.secs['general']['nodes'])
            if self.secs['general']['nodes'] < 0:
                raise ValueError
        except ValueError:
            fail('The values of the nodes is incorrect.')

        self.secs['general']['opt_mpi'] = True if self.secs['general']['opt_mpi'].lower() in ('yes', '1', 'true') else False

        try:
            self.secs['general']['maxwarn'] = int(self.secs['general']['maxwarn'])
            if self.secs['general']['maxwarn'] < 0:
                raise ValueError
        except ValueError:
            fail('The values of the maxwarn is incorrect.')

        self.secs['general']['email'] = self.secs['general'].get('email', '')
        self.secs['general']['title'] = self.secs['general'].get('title', os.path.abspath(os.path.curdir)[1:].replace('/', '_'))

    def __init_input(self):
        '''Initialize keys in [input]'''
        try:
            self.secs['input']['box'] = tuple(float(i) for i in self.secs['input']['box'].split())
        except ValueError:
            fail('The values of box size are incorrect. They are must be float.')
        if len(self.secs['input']['box']) != 3:
            fail('The number of the box size value must be 3.')

    def __init_component(self):
        '''Initialize keys in [component]'''
        self.secs['component']['itps'] = tuple(self.secs['component']['itps'].split())
        try:
            self.secs['component']['mols_num'] = tuple(int(i) for i in self.secs['component']['mols_num'].split())
            self.secs['component']['mols_name'] = tuple(i for i in self.secs['component']['mols_name'].split())
        except ValueError:
            fail('The values of the number of molecule(mols_num) are incorrect. They are must be int.')
        if len(self.secs['component']['itps']) != len(self.secs['component']['mols_num']) != len(self.secs['component']['mols_name']):
            fail('The number of itps, mols_num and mols_name in [component] is not equal.')

        self.secs['component']['residures'] = tuple(self.secs['component']['residures'].split())
        if len(self.secs['component']['residures']) <= 0:
            fail('The number of residures in system should be bigger than 0')
        self.secs['component']['extra_itps'] = tuple(self.secs['component']['extra_itps'].split())

    def __init_ions(self):
        '''Initialize keys in [component]'''
        self.secs['ions']['name'] = tuple(self.secs['ions']['name'].upper().split())
        try:
            self.secs['ions']['num'] = tuple(int(i) for i in self.secs['ions']['num'].split())
        except ValueError:
            fail('The values of ions number are incorrect. They are must be int.')

        if len(self.secs['ions']['name']) != len(self.secs['ions']['num']):
            fail('The number of values of ions name and ions num in [ions] is not equal.')

    def __init_pr(self):
        '''Initialize keys in [pr]'''
        self.secs['pr']['residures'] = tuple(self.secs['pr']['residures'].split())
        self.secs['pr']['itps'] = tuple(self.secs['pr']['itps'].split())
        self.secs['pr']['pdbs'] = tuple(self.secs['pr']['pdbs'].split())
        if len(self.secs['pr']['residures']) != len(self.secs['pr']['itps']) != len(self.secs['pr']['pdbs']):
            fail('The number of residures, itps and pdbs in [pr] is not equal.')

    def __init_md(self):
        '''Initialize keys in [md]'''
        self.secs['md']['ems'] = tuple(self.secs['md']['ems'].split())
        self.secs['md']['mds'] = tuple(self.secs['md']['mds'].split())
        if len(self.secs['md']['ems']) == 0:
            fail('You must supply at least one step of energy minimization in ems option in [md] section.')
        if len(self.secs['md']['mds']) == 0:
            fail('You must supply at least one step of md simulation in mds option in [md] section.')

    def __init_analysis(self):
        '''Initialize options in [analysis]

        The rdf will be parsered as:
            [[{id:type,...}, {id:type, ...}], ...]
            id is a sign for group, such as HW1, HW2, OW, SDmso etc.
            type is t or a, t is type for group, a is name for group, e.g. HW1 is one of the names of atom H in spce, and another name is HW2, and the both type is HW.
            The first () continues groups for rdf analysis. The second () is the only two groups. The dict contains groups will be joined as one group in ndx file.
        '''
        # rdf
        rdf = []
        one = [i.split() for i in self.secs['analysis']['rdf'].split(',')]
        for i in one:
            inner_1 = []
            for j in i:
                two = j.split('|')
                inner_2 = OrderedDict()
                for k in two:
                    if k.startswith('a#'):
                        inner_2[k[2:]] = 'a'
                    else:
                        inner_2[k] = 't'
                inner_1.append(inner_2)
            if len(inner_1) > 0:
                rdf.append(tuple(inner_1))
        self.secs['analysis']['rdf'] = tuple(rdf)

    def __init_simulation(self, section):
        '''
        Initialize sections or simulations in [md] ems and mds values.

        This function will do these things:
            1. Check the value wheather it is reasonable.
            2. Copy the default mdp parameters from default section such as [em_default] and [md_default] to this section dict.
            3. Initialize the value in this section making it can be used directly.
        '''
        type = self.get(section, 'type', fallback='')
        if type in ('steep', 'cg'):
            self.ems[section] = {k: v for k, v in self.items('em_default')}
            self.ems[section].update({k: v for k, v in self.items(section)})
            self.__init_userdefine_option(section, self.ems[section])
        elif type in ('npt', 'nvt', 'npt_pr', 'nvt_pr'):
            self.mds[section] = {k: v for k, v in self.items('md_default')}
            self.mds[section].update({k: v for k, v in self.items(section)})
            if type[-3:] == '_pr':
                self.__init_mds_pr(section)
            if type[:3] == 'npt':
                self.__init_mds_npt(section)
            self.__init_common_mds(section)
            self.__init_userdefine_option(section, self.mds[section])
        elif type == 'ext':
            self.mds[section] = {k: v for k, v in self.items(section)}
            self.__init_mds_ext(section)
            self.__init_userdefine_option(section, self.mds[section])
        else:
            fail('The type in [{0}] section is incorrect.')

    def __init_userdefine_option(self, section, sec_dict):
        '''Initialize and check the user define option in ems or mds.'''
        # nodes
        try:
            sec_dict['nodes'] = int(sec_dict.get('nodes', 0))
        except ValueError:
            fail('The values of nodes in [{0}] section should be integer.')

        # mdp
        sec_dict['mdp'] = sec_dict.get('mdp', './include/' + section + '.mdp')
        # analysis_dir
        sec_dict['analysis_dir'] = sec_dict.get('analysis_dir', './result/' + section)

        # mail_results
        sec_dict['mail_results'] = sec_dict.get('mail_results', '').split()

    def __init_mds_pr(self, section):
        '''Initialize and check the position restraint option in mds.'''
        resi = self.mds[section].get('pr_resi')
        if not resi in self.secs['pr']['residures']:
            fail('You must supply a residure in residures options of [pr] section')
        self.mds[section]['define'] = '-DPOSRES_' + resi.upper()

    def __init_mds_npt(self, section):
        '''
        Initialize and check the npt keys&values in mds
        '''
        pcoupltype = self[section].get('pcoupltype')
        # the value_num of tau-p, ref-p, compressibility
        value_num = tuple(len(self[section].get(key, '').split()) for key in ('compressibility', 'ref-p'))
        if value_num.count(value_num[0]) != 2:
            fail('The number of values of tau-p, ref-p, compressibility should be equal in [{}] section.'.format(section))
        num = value_num[0]

        if not pcoupltype in ('semiisotropic', 'isotropic', 'anisotropic', 'surface-tension'):
            fail('The pcoupltype should be one of the semiisotropic, isotropic, anisotropic, surface-tension in [{}] section.'.format(section))
        elif pcoupltype.lower() == 'isotropic' and num != 1:
            fail('When pcoupltype is isotropic , there should be only 1 values of tau-p, ref-p, compressibility in [{}] section.'.format(section))
        elif (pcoupltype.lower() == 'semiisotropic' or pcoupltype.lower() == 'surface-tension') and num != 2:
            fail('When pcoupltype is semiisotropic or surface-tension, there should be 2 values of tau-p, ref-p, compressibility in [{}] section.'.format(section))
        elif pcoupltype.lower() == 'anisotropic' and num != 6:
            fail('When pcoupltype is semiisotropic, there should be 6 values of tau-p, ref-p, compressibility in [{}] section.'.format(section))

    def __init_mds_ext(self, sec):
        '''
        Initialize and check the ext options in mds
        '''
        mds = self.secs['md']['mds']
        if not sec in mds:
            fail('The ext simulation only can be used in mds not energy minimization in [{}] section.'.format(sec))
        for i in range(mds.index(sec) - 1, -1, -1):
            if self.mds[mds[i]]['type'] != 'ext':
                break
        else:
            fail('Before the ext simulation there must be a real simulation such as npt, nvt in [{}] section.'.format(sec))

        try:
            self.mds[sec]['time'] = self.getint(sec, 'time')
            if self.mds[sec]['time'] <= 0:
                raise ValueError
        except ValueError:
            fail('The values of time in [{0}] section should be integer and > 0.'.format(sec))

    def __init_common_mds(self, section):
        '''
        Initialize and check the common keys&values in ems and mds
        '''
        sec = self.mds[section]
        needed_option = ('nsteps', 'continuation', 'tcoupl', 'tau-t', 'ref-t', 'gen_vel')
        for option in needed_option:
            if not self.has_option(section, option):
                fail('There is no {0} option in [{1}] section.'.format(option, section))

        # check 'continuation' and 'gen_vel'
        for key in ('continuation', 'gen_vel'):
            if not sec[key] in ('yes', 'no'):
                fail('The value of {0} in [{1}] section should be yes or no'.format(key, section))
        if sec['continuation'] == sec['gen_vel']:
            fail('The value of continuation and gen_vel should not be {0} at the same time'.format(sec['continuation']))
        # Check the continuation when the section is the first md
        if self.secs['md']['mds'][0] == section and sec['continuation'] == 'yes':
            warning('Are you sure set continuation to be yes? The [{0}] section is the first md course.'.format(section))
        if self.secs['md']['mds'][0] != section and sec['continuation'] == 'no':
            warning('Are you sure set continuation to be no? The [{0}] section is the not the first md course.'.format(section))

        # check gen_temp
        if sec['gen_vel'] == 'yes':
            try:
                sec['gen_temp'] = float(sec.get('gen_temp', ''))
            except ValueError:
                fail('You must supply a proper value of gen_temp in [{0}] section.'.format(section))

        # handle tc_grps tau-t ref-t
        residures = list(self.secs['component']['residures'][:])
        if 'SOL' in residures:
            residures.remove('SOL')
            # If ions are added into system, the group should contain  water and all
            if len(self.secs['ions']['name']) > 0:
                # This group name need to be generated by make_ndx programme
                water_ions_group = 'SOL_' + '_'.join(self.secs['ions']['name'])
            else:
                water_ions_group = 'SOL'
            residures.append(water_ions_group)
        self.mds[section]['tc-grps'] = ' '.join(residures)
        self.mds[section]['tau-t'] = ((self.mds[section]['tau-t'] + ' ') * len(residures)).strip()
        self.mds[section]['ref-t'] = ((self.mds[section]['ref-t'] + ' ') * len(residures)).strip()
        if 'annealing_temp' in self.mds[section]:
            note('The tc-grps in {} section is {}. Please adjust your parameters properly. e.g. annealing_temp'.format(section, self.mds[section]['tc-grps']))

    def get_secitions_dict(self):
        '''Get the dict of the sections in ini file except the simulation section.'''
        return self.secs

    def get_ems_dict(self):
        '''Get the dict of energy minimization sections.'''
        return self.ems

    def get_mds_dict(self):
        '''Get the dict of molecular simulation sections.'''
        return self.mds


class CommandOut:
    '''Output the commands for running simulation'''
    ndx = 'system.ndx'
    gro = 'pdb.gro'

    def __init__(self, parser, output):
        self.output = output
        '''Get dicts from configuration parser'''
        self.secs = parser.get_secitions_dict()
        self.ems = parser.get_ems_dict()
        self.mds = parser.get_mds_dict()
        self.top = self.secs['input']['top'] + '.top'

    def __write(self, text, mode='a'):
        '''Write text to commands file.'''
        with open(self.output, mode) as f:
            f.write(text)

    def __cmd(self, cmd):
        '''Get the gromacs command according the precision.'''
        precision = self.secs['general']['precision']
        # adjust for new gromacs
        cmd = "gmx " + cmd
        if precision[:3] == 'mpi':
            cmd += '_mpi'
        if precision[-6:] == 'double':
            cmd += '_d'
        return cmd

    def generate(self, exec_analysis=False):
        '''输出命令'''
        # Create file
        self.__write('#!/bin/bash\n', 'w')

        # Convert start pdb
        self.__pdb_conv(self.secs['input']['start_pdb'] + '.pdb')

        # add ions
        if len(self.secs['ions']['name']) > 0:
            self.__genion()

        # Generate index file
        self.__make_ndx(self.gro)
        self.__genrestr()

        #生成能量最小化命令
        for sec in self.ems.keys():
            # grompp
            self.__grompp(sec)
            # mdrun
            self.__mdrun(sec)

        # 生成模拟过程的命令
        for i in self.mds.keys():
            md_type = self.mds[i]['type']
            if md_type == 'ext':
                # convert-tpr
                self.__convert_tpr(i)
            else:
                # grompp
                self.__grompp(i)
            # mdrun
            self.__mdrun(i)
            if len(self.mds[i]['mail_results']) > 0:
                self.__mail_energy_result(i)
            if exec_analysis:
                self.__exec_analysis(i)

    def __pdb_conv(self, pdb):
        '''Output the commands for converting pdb'''
        editconf = self.__cmd('editconf')
        box = ' '.join([str(i) for i in self.secs['input']['box']])
        self.__write('{0} -f {1} -box {2} -o {3}\n'.format(editconf, pdb, box, self.gro))

    def __genion(self):
        '''Output genion commands'''
        ions = self.secs['ions']['name']
        nums = self.secs['ions']['num']

        # Use the first em mdp file
        em_sec = self.secs['md']['ems'][0]
        mdp = self.ems[em_sec]['mdp']
        gro = 'pdb.gro'
        tpr = 'ions.tpr'
        grompp = self.__cmd('grompp')
        genion = self.__cmd('genion')
        positive_ions = ('NA', 'K', 'MG', 'CA', 'ZN')
        negative_ions = ('CL', 'BR', 'I', 'F')

        for i, ion in enumerate(ions):
            grompp_cmd = '{} -f {} -c {} -p {} -o {}'.format(grompp, mdp, gro, self.top, tpr)
            genion_cmd = 'echo SOL | {} -s {} -o {} -p {}'.format(genion, tpr, gro, self.top)
            if ion in positive_ions:
                genion_cmd += ' -pname ' + ion + ' -np ' + str(nums[i])
            elif ion in negative_ions:
                genion_cmd += ' -nname ' + ion + ' -nn ' + str(nums[i])
            else:
                fail('The ion name is incorrect. Now it only supports {}'.format(', '.join(list(positive_ions) + list(negative_ions))))
            self.__write(grompp_cmd + '\n' + genion_cmd + '\n')

        # neutralization
        grompp_cmd = '{} -f {} -c {} -p {} -o {}'.format(grompp, mdp, gro, self.top, tpr)
        genion_cmd = 'echo SOL | {} -s {} -o {} -p {}'.format(genion, tpr, gro, self.top)
        genion_cmd += ' -pname NA -nname CL -conc 0.00000000001 -neutral'
        self.__write(grompp_cmd + '\n' + genion_cmd + '\n')

    def __make_ndx(self, gro):
        '''Generate index file'''
        # Generate tpr file
        grompp = self.__cmd('grompp')
        mdp = self.ems[list(self.ems.keys())[0]]['mdp']
        tpr = 'make_ndx.tpr'
        grompp_cmd = '{} -f {} -c {} -p {} -o make_ndx.tpr'.format(grompp, mdp, gro, self.top)
        self.__write(grompp_cmd + '\n')

        ions = self.secs['ions']['name']
        syntax = '"'
        # Generate groups including SOL and ions
        if len(ions) > 0:
            syntax += 'r SOL|r ' + '|r '.join(ions) + r'\n'
        rdf = self.secs['analysis']['rdf']
        if len(rdf) > 0:
            groups = set(tuple(j.items()) for i in rdf for j in i)
            # remove ions, since there are also ions in ndx
            groups.difference_update(set(((i, 't'),) for i in ('NA', 'K', 'MG', 'CA', 'ZN', 'CL')))
            for group in tuple(groups):
                syntax += '|'.join([i[1] + ' ' + i[0] for i in group]) + r'\n'
        syntax += r'q\n"'
        self.__write('echo -e {}| {} -f {} -o {} \n'.format(syntax, self.__cmd('make_ndx'), tpr, self.ndx, tpr))

    def __genrestr(self):
        '''Output position restraint itp file'''
        genrestr = self.__cmd('genrestr')

        # The residures to be restraint
        resis = self.secs['pr']['residures']
        pdbs = self.secs['pr']['pdbs']

        for i, res in enumerate(resis):
            make_ndx_cmd = 'echo "q\\n"|{} -f {}.pdb -o {}.ndx'.format(self.__cmd('make_ndx'), pdbs[i], res.lower())
            # getnr.py is a python programm can get the index of a residure, it must be in the path.
            genrestr_cmd = 'echo `getnr.py {}.ndx {}`| {} -f {}.pdb -o posre_{}.itp -fc 1000 1000 1000'.format(res.lower(), res, genrestr, pdbs[i], res.lower())
            self.__write(make_ndx_cmd + '\n')
            self.__write(genrestr_cmd + '\n')

    def __convert_tpr(self, sec):
        '''Output convert-tpr command for ext'''
        tpr = sec + '.tpr'
        convert_tpr = self.__cmd('convert-tpr')
        time = self.mds[sec]['time']
        mds = self.secs['md']['mds']
        last_sec = mds[mds.index(sec) - 1] if mds.index(sec) > 0 else None
        last_tpr = last_sec + '.tpr'
        convert_tpr_cmd = '{} -s {} -extend {} -o {}'.format(convert_tpr, last_tpr, time, tpr)
        self.__write(convert_tpr_cmd + '\n')

    def __grompp(self, sec):
        '''Output grompp command.'''
        secs = self.ems if sec in self.secs['md']['ems'] else self.mds
        grompp = self.__cmd('grompp')
        mdp = secs[sec]['mdp']
        tpr = sec + '.tpr'

        # Previous section, if not exist, set it to be None
        keys = tuple(secs.keys())
        last_sec = keys[keys.index(sec) - 1] if keys.index(sec) > 0 else None
        if sec == list(secs.keys())[0]:
            gro = 'pdb.gro' if secs[sec]['type'] in ('steep', 'cg') else self.secs['md']['ems'][-1]
        else:
            gro = last_sec + '.gro'

        grompp_cmd = grompp + ' -f ' + mdp + ' -c ' + gro + ' -p ' + self.top + ' -o ' + tpr + ' -n ' + self.ndx

        # If previous md type is not energy minimization, there must be a -t option to get checkpoint from it
        if last_sec in self.secs['md']['mds']:
            grompp_cmd += ' -t ' + last_sec + '.cpt'

        # maxwarn option
        maxwarn = self.secs['general']['maxwarn']
        if maxwarn > 0:
            # If there is maxwarn, add -maxwarn option
            grompp_cmd += ' -maxwarn ' + str(maxwarn)

        # If opt_pme_load is off
        if sec in self.secs['md']['ems'] or not self.secs['md']['opt_pme_load'] in ['1', 'yes', 'Yes', 'YES']:
            self.__write(grompp_cmd + '\n')
        else:
            # It will use opt_pme.py programm
            note("It will use opt_pme.py programm, please insure the programm be in path.")
            grompp_cmd += ' 2>&1 | tee /tmp/mdtool_opt_pme.txt'
            get_load_value = 'pme_load=`cat /tmp/mdtool_opt_pme.txt|grep "Estimate for the relative computational load of the PME mesh part: "|awk -F":[[:space:]]" \'{print $2}\'`'
            pme_load_std = self.secs['md']['pme_load']
            circle = grompp_cmd + ';' + get_load_value + ';'
            circle += 'while [[ "$pme_load" != "%s" ]];' % pme_load_std
            # Optimize
            circle += 'do fourierspacing=`opt_pme.py $pme_load %s "%s"`' % (pme_load_std, mdp)
            # If fourierspacing is too big, as now is set to be 1.00, abandon.
            circle += ';if [[ ${fourierspacing:0:1} -eq 1 ]]; then break;fi'
            circle += ';rm mdout.mdp;rm %s.tpr;%s;%s' % (sec, grompp_cmd, get_load_value)
            circle += ';done\n'

            self.__write(circle)

    def __mdrun(self, sec):
        '''Output mdrun command'''
        #是否是用mpi
        is_mpi = True if self.secs['general']['precision'][:3] == 'mpi' else False
        if is_mpi:
            self.__mdrun_mpi(sec)
        else:
            self.__mdrun_single(sec)

    def __mdrun_single(self, sec):
        '''Single machine mdrun command.'''
        mdrun = self.__cmd('mdrun')
        mds = self.ems if sec in self.ems.keys() else self.mds
        md = mds[sec]
        # Previous section, if not exist, set it to be None
        keys = tuple(mds.keys())
        last_sec = keys[keys.index(sec) - 1] if keys.index(sec) > 0 else None

        mdrun_cmd = mdrun + ' -deffnm ' + sec
        nodes = self.secs['general']['nodes']
        if nodes > 0:
            mdrun_cmd += ' -nt ' + str(nodes)

        if md['type'] == 'ext' and last_sec:
            mdrun_cmd += ' -cpi ' + last_sec + '.cpt'
        # If nodes > 0, use -nt option
        nodes = md['nodes']
        if nodes > 0:
            mdrun_cmd += ' -nt ' + str(nodes)

        self.__write(mdrun_cmd + '\n')

    def __mdrun_mpi(self, sec):
        '''Output mpiexec mdrun command.'''
        mds = self.ems if sec in self.ems.keys() else self.mds
        md = mds[sec]
        keys = tuple(mds.keys())
        last_sec = keys[keys.index(sec) - 1] if keys.index(sec) > 0 else None

        mdrun = self.__cmd('mdrun')
        g_tune_pme = self.__cmd('tune_pme')
        tpr = sec + '.tpr'

        nodes = self.secs['general']['nodes']
        if nodes > 0:
            if sec in self.secs['md']['ems'] or self.secs['general']['opt_mpi'] == 0:
                md_cmd = 'mpiexec -np ' + str(nodes) + ' ' + mdrun
            else:
                md_cmd = g_tune_pme + ' -np ' + str(nodes) + ' -s ' + tpr + ' -launch'
        else:
            fail('The mpi nodes must be bigger than 0.')

        md_cmd += ' -deffnm ' + sec

        if md['type'] == 'ext':
            md_cmd += ' -cpi ' + last_sec + '.cpt'
        self.__write(md_cmd + '\n')

    def __mail_energy_result(self, sec):
        '''mail md energy analysising result to email which is in [general] section.'''
        count = len(self.mds[sec]['mail_results'])
        emailsNum = len(self.secs['general']['email'])
        if count > 0 and emailsNum > 0:
            for index, ener in enumerate(self.mds[sec]['mail_results']):
                self.__write('result{0}=`echo "{1}"| {2} -f {3}.edr -o /tmp/energy.xvg|grep "{4}"`;'.format(index, ener, self.__cmd('energy'), sec, ener.replace('*', r'\*')))
            self.__write(r'echo -e "I am `hostname`. Of the simulation "{}" the section {} is just finished at `date +"%Y-%m-%d %H:%M"`. The energy results are as following: \\n'.format(self.secs['general']['title'], sec))
            for i in range(count):
                self.__write(r'${result%s}\\n' % i)
            if len(self.secs['general']['email']) > 0:
                self.__write('"| mutt -s "Energy Result" %s > /dev/null;echo "Finished!"\n' % self.secs['general']['email'])

    def __exec_analysis(self, sec):
        '''Execute the analysis script after md process.'''
        md = self.mds[sec]
        script = md['analysis_dir'] + '_ass'
        self.__write('sh %s;echo "Finished!"\n' % script)


class AnalysisOut:
    '''Output analysis commands file for every simulation step'''
    ndx = 'system.ndx'

    def __init__(self, parser):
        '''Get dicts from configuration parser'''
        self.secs = parser.get_secitions_dict()
        self.ems = parser.get_ems_dict()
        self.mds = parser.get_mds_dict()

    def __write(self, filename, text, mode='a'):
        '''Write text to commands file.'''
        with open(filename, mode) as f:
            f.write(text)

    def __cmd(self, cmd):
        '''Get the gromacs command according the precision.'''
        precision = self.secs['general']['precision']
        if precision[:3] == 'mpi':
            cmd += '_mpi'
        if precision[-6:] == 'double':
            cmd += '_d'
        return cmd

    def output(self):
        '''output the analysis commands to file.'''
        mds = self.mds
        for sec in mds.keys():
            self.__common(sec)
            #if mds[sec]['type'][:3] == 'npt':
            #    self.__npt()

    def __common(self, sec):
        '''Output the common analysis commands for npt, nvt etc.'''
        mds = self.ems if sec in self.ems.keys() else self.mds
        md = mds[sec]

        # Previous section which is not ext
        md_list = list(mds.keys())
        for i in range(md_list.index(sec) - 1, -1, -1):
            if mds[md_list[i]]['type'] != 'ext':
                last_noext = mds[md_list[i]]
                break

        path = md['analysis_dir']
        # log
        log = path + '/analysis.log'
        ass = sec + '_ass'

        effect_md = md if md['type'] != 'ext' else last_noext
        try:
            if 'nstxtcout' in effect_md.keys() and int(effect_md['nstxtcout']) > 0:
                trr = 'xtc'
            else:
                trr = 'trr'
        except ValueError as err:
            fail('{}. The nstxtcout in [{}] section must be integer.'.format(err, sec))

        # Gromacs programm
        g_energy, g_rms, g_density, g_rdf, trjconv = (self.__cmd(i) for i in ("gmx energy", "gmx rms", "gmx density", "gmx rdf", "gmx trjconv"))

        # file name for analysis
        gro, edr, tpr, trr = (sec + i for i in (".gro", ".edr", ".tpr", ".trr"))

        # Begin outputing
        self.__write(ass, '#!/bin/bash', 'w')
        # Find the root directory of simulation
        self.__write(ass, '''
if [[ -f "../{0}" ]]
then
    cd ..
elif [[ -f "./{0}" ]]
then
    cd .
else
    echo -e "\\033[91mCan't find {0} file in current and parent directory. Please put this analysis command file to the root directory or its first subdirectory.\\033[0m"
    exit 1
fi'''.format(self.ndx))

        # Make directory
        self.__write(ass, '\n[[ -d "{0}" ]] || mkdir -pv {0}'.format(path))

        # Header comments
        self.__write(ass, '\necho --------------------{}-------------------- > {}\n'.format(ass, log))

        # Groups whose nrs are to be gotten
        rdf = self.secs['analysis']['rdf']
        residures = self.secs['component']['residures']
        groups = []
        if len(rdf) > 0:
            groups += ['_'.join(j.keys()).upper() for i in rdf for j in i]
        groups = list(set(groups))
        groups += ["System"] + list(residures) + list(self.secs['ions']['name'])
        for_groups = ' '.join(groups)
        get_nr = '''for i in {}
do
    eval $i=`getnr.py {} $i `
    if (( ${}!i{} < 0 ))
    then
        echo "The nr of group $i <0. Please check the residures and rdf in configuration file."
        exit 1
    fi
done'''.format(for_groups, self.ndx, '{', '}')
        self.__write(ass, get_nr + '\n')

        # Energy
        for energy in ['potential', 'total-energy', 'pressure', 'temperature', '#Surf*SurfTen', 'density', 'Pres-XX', 'Pres-YY', 'Pres-ZZ', 'Box-X', 'Box-Y', 'Box-Z']:
            self.__write(ass, 'echo "{0}"| {1} -f {2} -o {3}/{0}.xvg >> {4}\n'.format(energy, g_energy, edr, path, log))

        # density
        self.__write(ass, 'echo 0 | {0} -f {1} -n {2} -s {3} -d z -dens mass -o {4}/density_system.xvg >> {5} \n'.format(g_density, trr, self.ndx, tpr, path, log))
        meta = list(residures) + list(self.secs['ions']['name'])
        for density in meta:
            cmd = 'echo ${0} | {1} -f {2} -n {3} -s {4} -d z -dens mass -o {5}/density_{0}.xvg >> {6}\n'.format(density, g_density, trr, self.ndx, tpr, path, log)
            self.__write(ass, cmd)
        #rms
        for i in residures:
            self.__write(ass, 'echo ${0} ${0}|{1} -f {2} -s {3} -o {4}/rmsd_{0}.xvg >> {5}\n'.format(i, g_rms, trr, tpr, path, log))

        # rdf
        if len(rdf) > 0:
            rdf_group = [('_'.join(i[0].keys()), '_'.join(i[1].keys())) for i in rdf]
            for (i, j) in rdf_group:
                self.__write(ass, 'echo -e "${0}\\n${1}\\n" |{2} -f {3} -n {4} -o {5}/rdf_{0}_{1}.xvg >> {6}\n'.format(i.upper(), j.upper(), g_rdf, trr, self.ndx, path, log))


class TopOut:
    '''Output topology file.'''
    def __init__(self, parser):
        '''Get dicts from configuration parser'''
        self.secs = parser.get_secitions_dict()
        self.ems = parser.get_ems_dict()
        self.mds = parser.get_mds_dict()
        self.filename = self.secs['input']['top'] + '.top'

    def output(self):
        '''Output topology file'''
        ff = self.secs['input']['forcefield']
        itps = self.secs['component']['itps']
        #residures = self.secs['component']['residures']
        pr_itps = self.secs['pr']['itps']
        pr_residures = self.secs['pr']['residures']
        with open(self.filename, 'w') as f:
            f.write('#include "{0}.ff/forcefield.itp"\n'.format(ff))
            for i in itps:
                if i in ('spc', 'spce', 'tip3p', 'tip4p', 'tip5p'):
                    f.write('#include "{0}.ff/{1}.itp"\n'.format(ff, i))
                else:
                    f.write('#include "{}.itp"\n'.format(i))
                if i in pr_itps:
                    resi = pr_residures[pr_itps.index(i)]
                    # Position restraint
                    f.write('\n; Include Position restraint file\n#ifdef POSRES_{0}\n#include "posre_{1}.itp"\n#endif\n\n'.format(resi.upper(), resi.lower()))
            for j in self.secs['component']['extra_itps']:
                if j in ('ions',):
                    f.write('#include "{0}.ff/{1}.itp"\n'.format(ff, j))
                else:
                    f.write('#include "{}.itp"\n'.format(j))

            # [ system ]
            f.write('\n[ system ]\n; name\nsurfactants in water\n\n[ molecules ]\n; Compound\t#mols\n')
            for res, resnum in zip(self.secs['component']['mols_name'], self.secs['component']['mols_num']):
                f.write(res + '\t\t' + str(resnum) + '\n')


class MdpOut:
    '''Output the mdp file.'''
    def __init__(self, parser):
        '''Get dicts from configuration parser'''
        self.secs = parser.get_secitions_dict()
        self.ems = parser.get_ems_dict()
        self.mds = parser.get_mds_dict()

        # The option below is user define parameters, ignoring them
        self.ignore = ['type', 'nodes', 'mdp', 'analysis_dir', 'pr_resi', 'mail_results']

    def output(self):
        for i in self.ems.keys():
            self.__mdp(i)
        for i in self.mds.keys():
            if self.mds[i]['type'] != 'ext':
                self.__mdp(i)

    def __mdp(self, sec):
        md = self.ems[sec] if sec in self.ems.keys() else self.mds[sec]
        with open(sec + '.mdp', 'w') as f:
            for i in sorted(md.keys()):
                if i in self.ignore:
                    continue
                f.write('{}\t=\t{}\n'.format(i, md[i]))


def main():
    arg_parser = argparse.ArgumentParser(description='Read the configuration file and generate command file, analysis command file, topology file and mdp file.')
    arg_parser.add_argument('-c', '--command', action='store', help='Generate command file. Default file name is commands.)')
    arg_parser.add_argument('-a', '--analysis', action='store_true', help='Generate analysis command file.')
    arg_parser.add_argument('-t', '--top', action='store_true', help='Generate topology file.')
    arg_parser.add_argument('-m', '--mdp', action='store_true', help='Generate mdp files.')
    arg_parser.add_argument('-i', '--input', action='store', required=True, help='Configuration file. INI file is recommended.')
    arg_parser.add_argument('--exec-analysis', action='store_true', help='Execute analysis script after every md process.')
    args = arg_parser.parse_args()
    conf_f = open(args.input)
    reader = ConfigReader(conf_f)
    if args.command:
        command_ge = CommandOut(reader, args.command)
        command_ge.generate(args.exec_analysis)
    if args.analysis:
        analysis = AnalysisOut(reader)
        analysis.output()
    if args.top:
        top = TopOut(reader)
        top.output()
    if args.mdp:
        mdp = MdpOut(reader)
        mdp.output()


def warning(string):
    '''Give a colorful display for string as a warning'''
    warning = '\033[93m'
    endc = '\033[0m'
    print(warning + "Warning: " + string + endc)


def fail(string):
    '''Give a colorful display for string as a fail'''
    fail = '\033[91m'
    endc = '\033[0m'
    print(fail + "Error: " + string + endc)
    exit(1)


def note(string):
    '''Give a colorful display for string as a note'''
    fail = '\033[93m'
    endc = '\033[0m'
    print(fail + "Note: " + string + endc)

if __name__ == '__main__':
    main()
