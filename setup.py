#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
""" setup script to install spacepy

Authors
-------
The SpacePy Team
Los Alamos National Laboratory

Copyright ©2010 Los Alamos National Security, LLC.
"""

import os, sys, shutil, getopt, glob, re
import subprocess
from distutils.core import setup
from distutils.command.build import build as _build
from distutils.command.install import install as _install
from distutils.command.sdist import sdist as _sdist
import distutils.ccompiler
import distutils.dep_util
from distutils.dist import Distribution as _Distribution
import distutils.sysconfig
from distutils.errors import DistutilsOptionError
from os import environ as ENVIRON


def subst(pattern, replacement, filestr,
          pattern_matching_modifiers=None):
    """
    replace pattern by replacement in string
    pattern_matching_modifiers: re.DOTALL, re.MULTILINE, etc.
    """
    if pattern_matching_modifiers is not None:
        cp = re.compile(pattern, pattern_matching_modifiers)
    else:
        cp = re.compile(pattern)
    if cp.search(filestr):  # any occurrence of pattern?
        filestr = cp.sub(replacement, filestr)
    return filestr


def default_f2py():
    """Looks for f2py based on name of python executable
    Assumes any suffix to python should also apply to f2py.
    This picks up .exe, version numbers, etc.
    """
    interp = os.path.basename(sys.executable)
    if interp[0:6] == 'python':
        candidate = 'f2py' + interp[6:]
        for dir in os.environ['PATH'].split(os.pathsep):
            if os.path.isfile(os.path.join(dir, candidate)):
                return candidate
    if sys.platform == 'win32':
        return 'f2py.py'
    else:
        return 'f2py'


class build(_build):
    """Extends base distutils build to make pybats, libspacepy, irbem"""

    user_options = _build.user_options + [
        ('fcompiler=', None,
         'specify the fortran compiler to use: pg, gnu95, gnu [gnu95]'),
        ('f2py=', None,
         'specify name (or full path) of f2py executable [{0}]'.format(
        default_f2py())),
        ('build-docs', None,
         'Build documentation with Sphinx (default: copy pre-built) [False]'),
        ]

    def initialize_options(self):
        self.fcompiler = None
        self.f2py = None
        self.build_docs = None
        _build.initialize_options(self)

    def finalize_options(self):
        _build.finalize_options(self)
        install_obj = self.distribution.get_command_obj('install')
        if self.fcompiler == None:
            self.fcompiler = install_obj.fcompiler
            if self.fcompiler == None:
                self.fcompiler = 'gnu95'
        if not self.fcompiler in ('pg', 'gnu', 'gnu95'):
            raise DistutilsOptionError(
                '--fcompiler must be pg, gnu, gnu95')
        if self.build_docs == None:
            self.build_docs = install_obj.build_docs
            if self.build_docs == None:
                self.build_docs = False
        if self.f2py == None:
            self.f2py = install_obj.f2py
            if self.f2py == None:
                self.f2py = default_f2py()
        if self.compiler == None:
            self.compiler = install_obj.compiler
            if self.compiler == None:
                if sys.platform == 'win32':
                    self.compiler = 'mingw32'

    def compile_pybats(self):
        outdir = os.path.join(self.build_lib, 'spacepy', 'pybats')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outpath = os.path.join(outdir, 'ctrace2d.so')
        srcdir = os.path.join('spacepy', 'pybats')
        srcpaths = [os.path.join(srcdir, f)
                    for f in ('ctrace2d.pyf', 'trace2d.c')]
        if distutils.dep_util.newer_group(srcpaths, outpath):
            os.chdir(srcdir)
            try:
                os.system(
                    '{0} -c ctrace2d.pyf trace2d.c'.format(
                    self.f2py))
                outpath = os.path.join('..', '..', outpath)
                if os.path.exists(outpath):
                    os.remove(outpath)
                shutil.move('ctrace2d.so',
                            os.path.join('..', '..', outdir))
            except:
                self.distribution.add_warning(
                    'pybats compile failed; pybats will not be available.')
            finally:
                os.chdir('..')
                os.chdir('..')

    def compile_LANLstar(self):
        outdir = os.path.join(self.build_lib, 'spacepy', 'LANLstar')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outpath = os.path.join(outdir, 'libLANLstar.so')
        srcdir = os.path.join('spacepy', 'LANLstar')
        srcpath = os.path.join(srcdir, 'LANLstar.f')
        if distutils.dep_util.newer(srcpath, outpath):
            os.chdir(srcdir)
            try:
                os.system(
                    '{0} -c  LANLstar.f -m libLANLstar --fcompiler={1}'.format(
                    self.f2py, self.fcompiler))
                outpath = os.path.join('..', '..', outpath)
                if os.path.exists(outpath):
                    os.remove(outpath)
                shutil.move('libLANLstar.so',
                            os.path.join('..', '..', outdir))
            except:
                self.distribution.add_warning(
                    'LANLstar compile failed; LANLstar will not be available.')
            finally:
                os.chdir('..')
                os.chdir('..')

    def compile_irbempy(self):
        # 64 bit or 32 bit?"
        bit = len('%x' % sys.maxsize)*4
        fcompiler = self.fcompiler
        irbemdir = 'irbem-lib-2010-12-21-rev275'
        srcdir = os.path.join('spacepy', 'irbempy', irbemdir, 'source')
        outdir = os.path.join(os.path.abspath(self.build_lib),
                              'spacepy', 'irbempy')
        if sys.platform == 'win32':
            libfile = 'irbempylib.pyd'
        else:
            libfile = 'irbempylib.so'
        sofile = os.path.join(outdir, libfile)
        sources = glob.glob(os.path.join(srcdir, '*.f')) + \
                  glob.glob(os.path.join(srcdir, '*.inc'))
        if not distutils.dep_util.newer_group(sources, sofile):
            #up to date
            return
        if not sys.platform in ('darwin', 'linux2', 'win32'):
            self.distribution.add_warning(
                '%s not supported at this time. ' % sys.platform + 
                'IRBEM will not be available')
            return
        if self.fcompiler == 'pg' and sys.platform == 'darwin':
            self.distribution.add_warning(
                'Portland Group compiler "pg" not supported on Mac OS\n'
                'IRBEM will not be available.')
            return
        if self.fcompiler != 'gnu95' and sys.platform == 'win32':
            self.distribution.add_warning(
                'Only supported compiler on Win32 is gnu95\n'
                'IRBEM will not be available.')
            return

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        builddir = os.path.join(self.build_temp, 'irbem')
        if os.path.exists(builddir):
            shutil.rmtree(builddir)
        shutil.copytree(os.path.join('spacepy', 'irbempy', irbemdir),
                        builddir)
        # compile irbemlib
        olddir = os.getcwd()
        os.chdir(builddir)
        F90files = ['source/onera_desp_lib.f', 'source/CoordTrans.f', 'source/AE8_AP8.f']
        functions = ['make_lstar1', 'make_lstar_shell_splitting1', \
                     'coord_trans1','find_magequator1', 'find_mirror_point1',
                     'get_field1', 'get_ae8_ap8_flux', 'fly_in_nasa_aeap1',
                     'trace_field_line2_1', 'trace_field_line_towards_earth1']

        # call f2py
        os.system('{0} --overwrite-signature -m irbempylib -h irbempylib.pyf '
                  '{1} only: {2} :'.format(
            self.f2py, ' '.join(F90files), ' '.join(functions)))
        # intent(out) substitute list
        outlist = ['lm', 'lstar', 'blocal', 'bmin', 'xj', 'mlt', 'xout', 'bmin', 'posit', \
                   'xgeo', 'bmir', 'bl', 'bxgeo', 'flux', 'ind']
        inlist = ['sysaxesin', 'sysaxesout', 'iyr', 'idoy', 'secs', 'xin']
        fln = 'irbempylib.pyf'
        print('Substituting fortran intent(in/out) statements')
        with open(fln, 'r') as f:
            filestr = f.read()
        for item in inlist:
            filestr = subst( ':: '+item, ', intent(in) :: '+item, filestr)
        for item in outlist:
            filestr = subst( ':: '+item, ', intent(out) :: '+item, filestr)
        with open(fln, 'w') as f:
            f.write(filestr)

        # compile (platform dependent)
        os.chdir('source')
        compile_cmd32 = {
            'pg': 'pgf77 -c -Mnosecond_underscore -w -fastsse -fPIC *.f',
            'gnu': 'g77 -c -w -O2 -fPIC -fno-second-underscore *.f',
            'gnu95': 'gfortran -m32 -c -w -O2 -fPIC -ffixed-line-length-none *.f',
            }
        compile_cmd64 = {
            'pg': 'pgf77 -c -Mnosecond_underscore -w -fastsse -fPIC *.f',
            'gnu': 'g77 -c -w -O2 -fPIC -fno-second-underscore *.f',
            'gnu95': 'gfortran -m64 -c -w -O2 -fPIC -ffixed-line-length-none *.f',
            }
        f2py_flags = '--fcompiler={0}'.format(fcompiler)
        if fcompiler == 'gnu':
            f2py_flags += ' --f77flags=-fno-second-underscore'
        if self.compiler:
            f2py_flags += ' --compiler={0}'.format(self.compiler)
        if bit == 32:
            os.system(compile_cmd32[fcompiler])
        else:
            os.system(compile_cmd64[fcompiler])
        if sys.platform == 'darwin':
            os.system('libtool -static -o libBL2.a *.o')
        elif sys.platform == 'linux2':
            os.system('ar -r libBL2.a *.o')
            os.system('ranlib libBL2.a')
        elif sys.platform == 'win32':
            os.system('ar -r libBL2.a *.o')
            os.system('ranlib libBL2.a')            
        os.chdir('..')
        os.system(
            '{0} -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 '
            '{1}'.format(
            self.f2py, f2py_flags))
        try:
            shutil.move(libfile, sofile)
        except:
            self.distribution.add_warning(
                'irbemlib compile failed. '
                'Try a different Fortran compiler? (--fcompiler)')
        os.chdir(olddir)
        return

    def compile_libspacepy(self):
        """Compile the C library, libspacepy. JTN 20110224"""
        srcdir = os.path.join('spacepy', 'libspacepy')
        outdir = os.path.join(self.build_lib, 'spacepy')
        try:
            comp = distutils.ccompiler.new_compiler(compiler=self.compiler)
            if hasattr(distutils.ccompiler, 'customize_compiler'):
                distutils.ccompiler.customize_compiler(comp)
            else:
                distutils.sysconfig.customize_compiler(comp)
            sources = list(glob.glob(os.path.join(srcdir, '*.c')))
            objects = [os.path.join(self.build_temp, s[:-2] + '.o')
                       for s in sources]
            headers = list(glob.glob(os.path.join(srcdir, '*.h')))
            #Assume every .o file associated with similarly-named .c file,
            #and EVERY header file
            outdated = [s for s, o in zip(sources, objects)
                        if distutils.dep_util.newer(s, o) or
                        distutils.dep_util.newer_group(headers, o)]
            if outdated:
                comp.compile(outdated, output_dir=self.build_temp)
            libpath = os.path.join(
                outdir, comp.library_filename('spacepy', lib_type='shared'))
            if distutils.dep_util.newer_group(objects, libpath):
                comp.link_shared_lib(objects, 'spacepy', libraries=['m'],
                                     output_dir=outdir)
        except:
            self.distribution.add_warning(
                'libspacepy compile failed; some operations may be slow.')
            print('libspacepy compile failed:')
            (t, v, tb) = sys.exc_info()
            print(v)

    def copy_docs(self):
        """Copy documentation from pre-build Doc directory."""
        outdir = os.path.join(os.path.abspath(self.build_lib),
                              'spacepy', 'Doc')
        indir = os.path.join('Doc', 'build', 'html')
        if os.path.exists(outdir):
            return
        if not os.path.exists(indir):
            print("No pre-built documentation, attempting to build...")
            self.make_docs()
            return
        shutil.copytree(indir, outdir)

    def make_docs(self):
        """Create/update documentation with Sphinx."""
        try:
            import sphinx
            import numpydoc
        except:
            if self.build_docs:
                self.distribution.add_warning(
                "Numpydoc and sphinx required to build documentation.\n"
                "Help will not be available; try without --build-docs.")
            else:
                self.distribution.add_warning(
                "Numpydoc and sphinx required to build documentation.\n"
                "Help will not be available.")
        builddir = os.path.join(os.path.join(self.build_temp, 'doctrees'))
        indir = os.path.join('Doc', 'source')
        outdir = os.path.join(os.path.abspath(self.build_lib),
                              'spacepy', 'Doc')
        cmd = 'sphinx-build -b html -d {0} {1} {2}'.format(
            builddir, indir, outdir)
        subprocess.check_call(cmd.split())

    def run(self):
        """Actually perform the build"""
        self.compile_irbempy()
        self.compile_LANLstar()
        self.compile_pybats()
        self.compile_libspacepy()
        if self.build_docs:
            self.make_docs()
        else:
            self.copy_docs()
        _build.run(self)


class install(_install):
    """Extends base distutils install to check versions, install .spacepy"""

    user_options = _install.user_options + [
        ('fcompiler=', None,
         'specify the fortran compiler to use: pg, gnu95, gnu [gnu95]'),
        ('f2py=', None,
         'specify name (or full path) of f2py executable [{0}]'.format(
        default_f2py())),
        ('build-docs', None,
         'Build documentation with Sphinx (default: copy pre-built) [False]'),
        ]

    def initialize_options(self):
        self.fcompiler = None
        self.f2py = None
        self.build_docs = False
        self.compiler = None
        _install.initialize_options(self)

    def finalize_options(self):
        if self.fcompiler == None:
            self.fcompiler = 'gnu95'
        if not self.fcompiler in ('pg', 'gnu', 'gnu95'):
            raise DistutilsOptionError(
                '--fcompiler must be pg, gnu, gnu95')
        if self.f2py == None:
            self.f2py = default_f2py()
        _install.finalize_options(self)
        if self.compiler == None:
            if sys.platform == 'win32':
                self.compiler = 'mingw32'

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
            self.distribution.add_warning(
                'SciPy and matplotlib were not found.'
                'They are required for large parts of SpacePy.')
        else:
            print('Dependencies OK.')

        _install.run(self)


class sdist(_sdist):
    """Rebuild the docs before making a source distribution"""

    def run(self):
        builddir = os.path.join(os.path.join('Doc', 'build', 'doctrees'))
        indir = os.path.join('Doc', 'source')
        outdir = os.path.join('Doc', 'build', 'html')
        cmd = 'sphinx-build -b html -d {0} {1} {2}'.format(
            builddir, indir, outdir)
        subprocess.check_call(cmd.split())
        os.chdir('Doc')
        try:
            cmd = 'make latexpdf'
            subprocess.check_call(cmd.split())
            shutil.move(os.path.join('build', 'latex', 'SpacePy.pdf'),
                        '.')
        except:
            self.distribution.add_warning(
                'PDF documentation rebuild failed.')
            print('PDF documentation rebuild failed:')
            (t, v, tb) = sys.exc_info()
            print(v)
        finally:
            os.chdir('..')
        _sdist.run(self)


class Distribution(_Distribution):
    """Subclass of main distutils Distribution that adds support for warnings"""

    def add_warning(self, msg):
        """Add a warning for this instance of setup"""
        self._warnings.append(msg)

    def print_warnings(self):
        """Print out warnings from this execution of setup"""
        if not self._warnings:
            return
        print('\nsetup produced the following warnings. '
              'Some functionality may be missing.\n')
        for w in self._warnings:
            print(w)

    def run_commands(self):
        """Run all setup commands"""
        self._warnings = []
        _Distribution.run_commands(self)
        self.print_warnings()


pkg_files = ['irbempy/*.py', 'LANLstar/*.py',
    'pybats/*.py', 'pybats/*.out', 'pycdf/*.py', 'data/*']


# run setup from distutil
setup(name='spacepy',
      version='0.1',
      description='SpacePy: Tools for Space Science Applications',
      author='Steve Morley, Josef Koller, Dan Welling, Brian Larsen, Mike Henderson, Jon Niehof',
      author_email='spacepy@lanl.gov',
      url='http://www.spacepy.lanl.gov',
      requires=['numpy','scipy','matplotlib (>=0.99)'],
      packages=['spacepy',],
      package_data={'spacepy': pkg_files},
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI-Approved Open Source :: Python Software Foundation License',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Physics',
          'Topic :: Software Development :: Libraries :: Python Modules'
          ],
      cmdclass={'build': build,
                'install': install,
                'sdist': sdist,
                },
      distclass=Distribution,
      )
