#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# setup.py to install spacepy

__version__ = "$Revision: 1.48 $, $Date: 2011/03/01 22:58:18 $"
__author__ = 'The SpacePy Team, Los Alamos National Lab (spacepy@lanl.gov)'

import os, sys, shutil, getopt, warnings
from distutils.core import setup
from distutils.command.build import build as _build
from distutils.command.install import install as _install
import distutils.ccompiler
import distutils.dep_util
from distutils.errors import DistutilsOptionError
import distutils.sysconfig
import glob
from os import environ as ENVIRON
import re

# -------------------------------------
def subst(pattern, replacement, filestr,
          pattern_matching_modifiers=None):
          
    """
    replace pattern by replacement in file
    pattern_matching_modifiers: re.DOTALL, re.MULTILINE, etc.
    """
    
    import re, shutil
    
    
    if pattern_matching_modifiers is not None:
        cp = re.compile(pattern, pattern_matching_modifiers)
    else:
        cp = re.compile(pattern)

    if cp.search(filestr):  # any occurence of pattern?
        filestr = cp.sub(replacement, filestr)
        
    return filestr

# -------------------------------------


#import tooblox by reading file from repository
# this will provide mostly the query_yes_no function and update()
exec(compile(open('spacepy/toolbox.py').read(), 'spacepy/toolbox.py', 'exec'))

class build(_build):
    """Extends base distutils build to make pybats, libspacepy, irbem"""

    user_options = _build.user_options + [
        ('fcompiler=', None,
         'specify the fortran compiler to use: pg, gnu95, gnu [gnu95]'),
        ]

    def initialize_options(self):
        self.fcompiler = None
        _build.initialize_options(self)

    def finalize_options(self):
        _build.finalize_options(self)
        if self.fcompiler == None:
            self.fcompiler = \
                 self.distribution.get_command_obj('install').fcompiler
            if self.fcompiler == None:
                self.fcompiler = 'gnu95'
        if not self.fcompiler in ('pg', 'gnu', 'gnu95'):
            raise DistutilsOptionError(
                '--fcompiler must be pg, gnu, gnu95')

    # -------------------------------------
    def compile_pybats(self):
        os.chdir(os.path.join('spacepy', 'pybats'))
        os.system('f2py -c ctrace2d.pyf trace2d.c')
        os.chdir('..')
        os.chdir('..')
    
    # -------------------------------------
    def compile_irbempy(self):
        if not sys.platform in ('darwin', 'linux2'):
            print('%s not supported at this time' % sys.platform)
            print('IRBEM will not be available')
            return 

        fcompiler = self.fcompiler
        # compile irbemlib
        os.chdir('spacepy/irbempy/irbem-lib-2010-12-21-rev275')

        F90files = ['source/onera_desp_lib.f', 'source/CoordTrans.f', 'source/AE8_AP8.f']
        functions = ['make_lstar1', 'make_lstar_shell_splitting1', \
                     'coord_trans1 find_magequator1', 'find_mirror_point1', 
                     'get_field1', 'get_ae8_ap8_flux', 'fly_in_nasa_aeap1', 
                     'trace_field_line2_1', 'trace_field_line_towards_earth1']

        # call f2py
        os.system('f2py --overwrite-signature -m irbempylib -h irbempylib.pyf '+' '.join(F90files) \
                  +' only: ' + ' '.join(functions) + ' :')
        # intent(out) substitute list
        outlist = ['lm', 'lstar', 'blocal', 'bmin', 'xj', 'mlt', 'xout', 'bmin', 'posit', \
                   'xgeo', 'bmir', 'bl', 'bxgeo', 'flux', 'ind']
        inlist = ['sysaxesin', 'sysaxesout', 'iyr', 'idoy', 'secs', 'xin']
        fln = 'irbempylib.pyf'
        print('Substituting fortran intent(in/out) statements')
        f = open(fln, 'r')
        filestr = f.read()
        f.close()
        for item in inlist:
            filestr = subst( ':: '+item, ', intent(in) :: '+item, filestr)
        for item in outlist:
            filestr = subst( ':: '+item, ', intent(out) :: '+item, filestr)
        f = open(fln, 'w')
        f.write(filestr)
        f.close()

        # compile (platform dependent)
        os.chdir('source')
        if sys.platform == 'darwin': # then mac OS
            if fcompiler == 'pg':
                raise NotImplementedError('Portland Group compiler option "pg" is not supported on Mac OS with f2py')
                #os.system('pgf77 -c -Mnosecond_underscore -w -fastsse *.f')
                #os.system('libtool -static -o libBL2.a *.o')
                #os.chdir('..')
                #os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=pg')
            elif fcompiler == 'gnu':
                os.system('g77 -c -w -O2 -fPIC -fno-second-underscore *.f')
                os.system('libtool -static -o libBL2.a *.o')
                os.chdir('..')
                os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu --f77flags=-fno-second-underscore ')
            else:
                os.system('gfortran -c -w -O2 -fPIC -ffixed-line-length-none *.f')
                os.system('libtool -static -o libBL2.a *.o')
                os.chdir('..')
                os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')
        elif sys.platform == 'linux2': # then linux
            if fcompiler == 'pg':
                os.system('pgf77 -c -Mnosecond_underscore -w -fastsse -fPIC *.f')
                os.system('ar -r libBL2.a *.o')
                os.system('ranlib libBL2.a')
                os.chdir('..')
                os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=pg')
            elif fcompiler == 'gnu':           
                os.system('g77 -c -w -O2 -fPIC -fno-second-underscore *.f')
                os.system('ar -r libBL2.a *.o')
                os.system('ranlib libBL2.a')
                os.chdir('..')
                os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu --f77flags=-fno-second-underscore')
            else:
                os.system('gfortran -c -w -O2 -fPIC -ffixed-line-length-none *.f')
                os.system('ar -r libBL2.a *.o')
                os.system('ranlib libBL2.a')
                os.chdir('..')
                os.system('f2py -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')

        err = os.system('mv -f irbempylib.so ../')
        if err:
            print '------------------------------------------------------'
            print 'WARNING: Something went wrong with compiling irbemlib.'
            print '------------------------------------------------------'
        os.chdir('../../..')
        return

    def compile_libspacepy(self):
        """Compile the C library, libspacepy. JTN 20110224"""
        olddir = os.getcwd()
        os.chdir(os.path.join('spacepy', 'libspacepy'))
        try:
            comp = distutils.ccompiler.new_compiler(compiler=self.compiler)
            extra = distutils.sysconfig.get_config_vars('CCSHARED')
            sources = list(glob.glob('*.c'))
            objects = [s[:-2] + '.o' for s in sources]
            headers = [s[:-2] + '.h' for s in sources
                       if os.path.exists(s[:-2] + '.h')]
            #Assume every .o file associated with similarly-named .c file,
            #and EVERY header file
            outdated = [s for s, o in zip(sources, objects)
                        if distutils.dep_util.newer(s, o) or
                        distutils.dep_util.newer_group(headers, o)]
            if outdated:
                comp.compile(outdated, extra_preargs=extra)
            if distutils.dep_util.newer_group(
                objects, comp.library_filename('spacepy', lib_type='shared')):
                comp.link_shared_lib(objects, 'spacepy', libraries=['m'])
        except:
            print('libspacepy compile failed; some operations may be slow:')
            (t, v, tb) = sys.exc_info()
            print(v)
        finally:
            os.chdir(olddir)

    def run(self):
        """Actually perform the build"""
        # run compile for irbem-lib first
        if os.path.exists(os.path.join('spacepy', 'irbempy', 'irbempylib.so')):
            ans = query_yes_no('\nDo you want to recompile the IRBEM library?', default="no")
            if ans=='yes':
                self.compile_irbempy()
        else:
            self.compile_irbempy()
        # Compile PyBats
        self.compile_pybats()
        # Compile libspacepy
        self.compile_libspacepy()
        _build.run(self)

class install(_install):
    """Extends base distutils install to check versions, install .spacepy"""

    user_options = _install.user_options + [
        ('fcompiler=', None,
         'specify the fortran compiler to use: pg, gnu95, gnu [gnu95]'),
        ]

    def initialize_options(self):
        self.fcompiler = None
        _install.initialize_options(self)

    def finalize_options(self):
        if self.fcompiler == None:
            self.fcompiler = 'gnu95'
        if not self.fcompiler in ('pg', 'gnu', 'gnu95'):
            raise DistutilsOptionError(
                '--fcompiler must be pg, gnu, gnu95')
        _install.finalize_options(self)

    def run(self):
        """Does all checks, perform install, makes .spacepy directory"""
        #test for python version 2.x where x>=6
        try:
            print('Checking Python >= 2.6...')
            dum = sys.version_info
            assert (dum[0] >= 2) or (dum[0] == 2 and dum[1] >= 6)
            print ('Checking for numpy...')
            import numpy
        except:
            raise Exception("""SpacePy requires Python 2.X, where X>=6.\n
            Numpy, Scipy and Matplotlib(>=0.99) are also required\n
            Please install suitable versions.""")
        try:
            print ('Checking for scipy...')
            import scipy
            print ('Checking for matplotlib...')
            import matplotlib
            assert matplotlib.compare_versions(
                matplotlib.__version__, '0.99.0')
        except:
            warnings.warn('''Missing Packages: SciPy and MatPlotLib are 
            required for large parts of this library.''')
        else:
            print('Dependencies OK.')

        _install.run(self)

        # create .spacepy in $HOME and move data
        # read-in .rc file first
        exec(compile(open('spacepy/data/spacepy.rc').read(), 'spacepy/data/spacepy.rc', 'exec'))
        if 'SPACEPY' in ENVIRON:
            DOT_FLN = os.path.join(ENVIRON['SPACEPY'], '.spacepy')
        else:
            DOT_FLN = os.path.join(ENVIRON['HOME'], '.spacepy')

        if os.path.exists(DOT_FLN):
            ans = query_yes_no('\n'+DOT_FLN+' already exists. Do you want to start fresh?', default="no")
            if ans=='no':
                fresh_install = False
            else:
                fresh_install = True
                i = 0
                while os.path.exists(DOT_FLN + '.bak.' + str(i)):
                    i += 1
                shutil.move(DOT_FLN, DOT_FLN+'.bak.'+str(i))
        else:
            fresh_install = True

        if fresh_install:
            os.mkdir(DOT_FLN)
            os.chmod(DOT_FLN, 0o777)
            datadir = os.path.join(DOT_FLN, 'data')
            os.mkdir(datadir)
            os.chmod(datadir, 0o777)
            shutil.copy(os.path.join('spacepy', 'data', 'spacepy.rc'), DOT_FLN)
            shutil.copy(os.path.join('spacepy', 'data', 'tai-utc.dat'),
                        datadir)

        # update/download packages
        if sys.version_info[0]<3:
            if fresh_install:
                dir = update()
                print("Data installed to " + dir)
            else:
                ans = query_yes_no("\nDo you want to update OMNI database and leap seconds table? (Internet connection required)", default = "no")
                if ans=='yes':
                    dir = update()
                else:
                    print("\nRemember to update OMNI and leap seconds table occasionally by running spacepy.toolbox.update()")
        else:
            print('''Updating OMNI and leap seconds on install is not currently supported for Python 3.X.''')
        print("\nThanks for installing SpacePy.")


pkg_files = ['irbempy/irbempylib.so', 'irbempy/*.py', 'doc/*.*', 'pybats/*.py', 'pybats/*.so'
    'pybats/*.out', 'pycdf/*.py', 'libspacepy/*spacepy*',]

# run setup from distutil
setup(name='spacepy',
      version='0.1',
      description='SpacePy: Tools for Space Science Applications',
      author='Steve Morley, Josef Koller, Dan Welling, Brian Larsen, Mike Henderson, Jon Niehof',
      author_email='spacepy@lanl.gov',
      url='http://www.spacepy.lanl.gov',
      requires=['numpy','scipy','matplotlib (>=0.99)'],
      packages=['spacepy','spacepy.sandbox'],
      package_data={'spacepy': pkg_files},
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Software Development :: Libraries :: Python Modules'
          ],
      cmdclass={'build': build,
                'install': install,
                },
      )
