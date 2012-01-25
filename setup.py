#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
""" setup script to install spacepy

Authors
-------
The SpacePy Team
Los Alamos National Laboratory

Copyright ©2010 - 2011 Los Alamos National Security, LLC.
"""

import os, sys, shutil, getopt, glob, re
import subprocess
from distutils.core import setup, Extension
from distutils.command.build import build as _build
from distutils.command.install import install as _install
from distutils.command.bdist_wininst import bdist_wininst as _bdist_wininst
from distutils.command.sdist import sdist as _sdist
import distutils.ccompiler
import distutils.dep_util
from distutils.dist import Distribution as _Distribution
import distutils.sysconfig
from distutils.errors import DistutilsOptionError

import numpy


#Patch out bad options in Python's view of mingw
if sys.platform == 'win32':
    import distutils.cygwinccompiler
    _Mingw32CCompiler = distutils.cygwinccompiler.Mingw32CCompiler
    class Mingw32CCompiler(_Mingw32CCompiler):
        def __init__(self, *args, **kwargs):
            _Mingw32CCompiler.__init__(self, *args, **kwargs)
            for executable in ('compiler', 'compiler_so', 'compiler_cxx',
                               'linker_exe', 'linker_so'):
                exe = getattr(self, executable)
                if '-mno-cygwin' in exe:
                    del exe[exe.index('-mno-cygwin')]
                    setattr(self, executable, exe)
    distutils.cygwinccompiler.Mingw32CCompiler = Mingw32CCompiler


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


def initialize_compiler_options(cmd):
    """Initialize the compiler options for a command"""
    cmd.fcompiler = None
    cmd.f2py = None
    cmd.compiler = None

    
def finalize_compiler_options(cmd):
    """Finalize compiler options for a command

    If compiler options (fcompiler, compiler, f2py) have not been
    specified for a command, check if they were specified for other
    commands on the command line--if so, grab from there. If not,
    set reasonable defaults.

    cmd: the command to finalize the options for
    """
    dist = cmd.distribution
    defaults = {'fcompiler': 'gnu95',
                'f2py': default_f2py(),
                'compiler': None}
    #Check all options on all other commands, reverting to default
    #as necessary
    for option in defaults:
        if getattr(cmd, option) == None:
            for c in dist.commands:
                other_cmd = dist.get_command_obj(c)
                if other_cmd == cmd or not hasattr(other_cmd, option):
                    continue
                if getattr(other_cmd, option) != None:
                    setattr(cmd, option, getattr(other_cmd, option))
                    break
            if getattr(cmd, option) == None:
                setattr(cmd, option, defaults[option])
    #Special-case defaults, checks
    if not cmd.fcompiler in ('pg', 'gnu', 'gnu95'):
        raise DistutilsOptionError(
            '--fcompiler must be pg, gnu, gnu95')
    if cmd.compiler == None and sys.platform == 'win32':
        cmd.compiler = 'mingw32'


compiler_options = [
        ('fcompiler=', None,
         'specify the fortran compiler to use: pg, gnu95, gnu [gnu95]'),
        ('f2py=', None,
         'specify name (or full path) of f2py executable [{0}]'.format(
        default_f2py())),
        ]


def rebuild_static_docs(dist=None):
    """Rebuild the 'static' documentation in Doc/build"""
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
    except:
        if dist != None:
            dist.add_warning(
                'PDF documentation rebuild failed.')
        print('PDF documentation rebuild failed:')
        (t, v, tb) = sys.exc_info()
        print(v)
    finally:
        os.chdir('..')


class build(_build):
    """Extends base distutils build to make pybats, libspacepy, irbem"""

    user_options = _build.user_options + compiler_options + [
        ('build-docs', None,
         'Build documentation with Sphinx (default: copy pre-built) [False]'),
        ]

    def initialize_options(self):
        self.build_docs = None
        initialize_compiler_options(self)
        _build.initialize_options(self)

    def finalize_options(self):
        _build.finalize_options(self)
        finalize_compiler_options(self)
        if self.build_docs == None:
            self.build_docs = self.distribution.get_command_obj(
                'install').build_docs
            if self.build_docs == None:
                self.build_docs = False

    def compile_LANLstar(self):
        outdir = os.path.join(self.build_lib, 'spacepy', 'LANLstar')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if sys.platform == 'win32':
            libfile = 'libLANLstar.pyd'
        else:
            libfile = 'libLANLstar.so'
        outpath = os.path.join(outdir, libfile)
        srcdir = os.path.join('spacepy', 'LANLstar')
        srcpath = os.path.join(srcdir, 'LANLstar.f')
        if distutils.dep_util.newer(srcpath, outpath):
            os.chdir(srcdir)
            try:
                f2py = self.f2py
                if self.compiler:
                    f2py += ' --compiler={0}'.format(self.compiler)
                os.system(
                    '{0} -c  LANLstar.f -m libLANLstar --fcompiler={1}'.format(
                    f2py, self.fcompiler))
                outpath = os.path.join('..', '..', outpath)
                if os.path.exists(outpath):
                    os.remove(outpath)
                shutil.move(libfile,
                            os.path.join('..', '..', outdir))
            except:
                self.distribution.add_warning(
                    'LANLstar compile failed; LANLstar will not be available.')
                print('LANLstar compile failed:')
                (t, v, tb) = sys.exc_info()
                print(v)
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
                return
            else:
                self.distribution.add_warning(
                "Numpydoc and sphinx required to build documentation.\n"
                "Help will not be available.")
                return
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
        self.compile_libspacepy()
        if self.build_docs:
            self.make_docs()
        else:
            self.copy_docs()
        _build.run(self)


class install(_install):
    """Extends base distutils install to check versions, install .spacepy"""

    user_options = _install.user_options + compiler_options + [
        ('build-docs', None,
         'Build documentation with Sphinx (default: copy pre-built) [False]'),
        ]

    def initialize_options(self):
        self.build_docs = False
        initialize_compiler_options(self)
        _install.initialize_options(self)

    def finalize_options(self):
        _install.finalize_options(self)
        finalize_compiler_options(self)

    def run(self):
        """Does all checks, perform install, makes .spacepy directory"""
        #test for python version 2.x where x>=6
        try:
            print('Checking Python >= 2.6...')
            dum = sys.version_info
            assert (dum[0] >= 2) or (dum[0] == 2 and dum[1] >= 6)
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
        #Files will be deleted in the order specified, so list files
        #before directories containing them!
        #Paths are relative to spacepy. Unix path separators are OK
        #Don't forget to delete the .pyc
        deletefiles = ['toolbox.py', 'toolbox.pyc']
        for f in deletefiles:
            path = os.path.join(self.install_lib, 'spacepy',
                                os.path.normpath(f)) #makes pathing portable
            if os.path.exists(path):
                print('Deleting {0} from old version of spacepy.'.format(path))
                if os.path.isdir(path):
                    os.rmdir(path)
                else:
                    os.remove(path)


class bdist_wininst(_bdist_wininst):
    """Handle compiler options for build on Windows install"""

    user_options = _bdist_wininst.user_options + compiler_options

    def initialize_options(self):
        initialize_compiler_options(self)
        _bdist_wininst.initialize_options(self)

    def finalize_options(self):
        _bdist_wininst.finalize_options(self)
        finalize_compiler_options(self)

    def copy_fortran_libs(self):
        """Copy the fortran runtime libraries into the build"""
        fortdir = None
        fortnames = None
        for p in os.environ['PATH'].split(';'):
            fortnames = [f for f in os.listdir(p)
                         if f[-4:].lower == '.dll' and
                         (f[:11] == 'libgfortran' or
                          f[:8] == 'libgcc_s' or
                          f[:11] == 'libquadmath')]
            if len(fortnames) == 3:
                fortdir = p
                break
        if fortdir is None:
            raise RuntimeError("Can't locate fortran libraries.")
        outdir = os.path.join(self.build_lib, 'spacepy', 'mingw')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        for f in fortnames:
            shutil.copy(os.path.join(fortdir, f), outdir)
            
    def run(self):
        self.copy_fortran_libs()
        rebuild_static_docs(self.distribution)
        _bdist_wininst.run(self)


class sdist(_sdist):
    """Rebuild the docs before making a source distribution"""

    def run(self):
        rebuild_static_docs(self.distribution)
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


packages = ['spacepy', 'spacepy.irbempy', 'spacepy.LANLstar',
            'spacepy.pycdf', 'spacepy.plot', 'spacepy.pybats', 'spacepy.time', 
            'spacepy.toolbox']
#If adding to package_data, also put in MANIFEST.in
package_data = ['data/*', 'pybats/sample_data/*']
pybats_ext = Extension('spacepy.pybats.ctrace2d',
                       sources=['spacepy/pybats/ctrace2dmodule.c'],
                       include_dirs=[numpy.get_include()])
dates_ext = Extension('spacepy.time._dates',
                                sources=['spacepy/time/_datesmodule.c'], 
                       include_dirs=[numpy.get_include()])
toolbox_ext = Extension('spacepy.toolbox._toolbox',
                                sources=['spacepy/toolbox/_toolboxmodule.c'], 
                       include_dirs=[numpy.get_include()])
                       
# run setup from distutil
setup(name='spacepy',
      version='0.1.1',
      description='SpacePy: Tools for Space Science Applications',
      ext_modules=[pybats_ext, dates_ext, toolbox_ext],
      author='Steve Morley, Josef Koller, Dan Welling, Brian Larsen, Mike Henderson, Jon Niehof',
      author_email='spacepy@lanl.gov',
      url='http://www.spacepy.lanl.gov',
      requires=['numpy','scipy','matplotlib (>=0.99)'],
      packages=packages,
      package_data={'spacepy': package_data},
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
                'bdist_wininst': bdist_wininst,
                'sdist': sdist,
                },
      distclass=Distribution,
      )
