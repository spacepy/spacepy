#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
""" setup script to install spacepy

Authors
-------
The SpacePy Team
Los Alamos National Laboratory

Copyright 2010 - 2013 Los Alamos National Security, LLC.
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


#These are files that are no longer in spacepy (or have been moved)
#having this here makes sure that during an upgrade old versions are
#not hanging out
#Files will be deleted in the order specified, so list files
#before directories containing them!
#Paths are relative to spacepy. Unix path separators are OK
#Don't forget to delete the .pyc
deletefiles = ['toolbox.py', 'toolbox.pyc', 'LANLstar/LANLstar.py',
               'LANLstar/LANLstar.pyc', 'LANLstar/libLANLstar.so',
               'LANLstar/LANLstar.pyd', 'LANLstar/__init__.py',
               'LANLstar/__init__.pyc', 'LANLstar',
               'time/__init__.py', 'time/__init__.pyc',
               'time/_dates.so', 'time/_dates.dylib',
               'time/_dates.pyd',
               'time/time.py', 'time/time.pyc', 'time',
               'data/LANLstar/*.net',
               'pycdf/_pycdf.*', 'toolbox/toolbox.py*']


def delete_old_files(basepath):
    """Delete files from old versions of spacepy, under a particular path"""
    for f in deletefiles:
        path = os.path.join(basepath, 'spacepy',
                            os.path.normpath(f)) #makes pathing portable
        for p in glob.glob(path):
            if os.path.exists(p):
                print('Deleting {0} from old version of spacepy.'.format(p))
                if os.path.isdir(p):
                    os.rmdir(p)
                else:
                    os.remove(p)


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
    interpdir, interp = os.path.split(sys.executable)
    if interp[0:6] == 'python':
        suffixes = [interp[6:], '-' + interp[6:]]
        if '.' in interp[6:]: #try slicing off suffix-of-suffix (e.g., exe)
            suffix = interp[6:-(interp[::-1].index('.') + 1)]
            suffixes.extend([suffix, '-' + suffix])
        candidates = ['f2py' + s for s in suffixes]
        for candidate in candidates:
            for c in [candidate, candidate + '.py']:
                for d in os.environ['PATH'].split(os.pathsep):
                    if os.path.isfile(os.path.join(d, c)):
                        return c
                    if os.path.isfile(os.path.join(interpdir, c)):
                        return os.path.join(interpdir, c) #need full path
    if sys.platform == 'win32':
        return 'f2py.py'
    else:
        return 'f2py'


def f2py_environment(fcompiler):
    """Prepare an OS environment for the f2py call
    This puts in the shared options if LDFLAGS is set
    """
    if not 'LDFLAGS' in os.environ:
        return None
    env = os.environ.copy()
    import numpy.distutils.fcompiler
    numpy.distutils.fcompiler.load_all_fcompiler_classes()
    if not fcompiler in numpy.distutils.fcompiler.fcompiler_class:
        return None #and hope for the best
    fcomp = numpy.distutils.fcompiler.fcompiler_class[fcompiler][1]()
    fcomp.customize()
    env['LDFLAGS'] = ' '.join(fcomp.get_flags_linker_so())
    return env


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
    if not cmd.fcompiler in ('pg', 'gnu', 'gnu95', 'intelem', 'intel', 'none'):
        raise DistutilsOptionError(
            '--fcompiler must be pg, gnu, gnu95, intelem, intel, none')
    if len('%x' % sys.maxsize)*4 == 32 and cmd.fcompiler == 'intelem':
        raise DistutilsOptionError(
            '--fcompiler=intelem requires a 64-bit architecture')
    if cmd.compiler == None and sys.platform == 'win32':
        cmd.compiler = 'mingw32'


compiler_options = [
        ('fcompiler=', None,
         'specify the fortran compiler to use: pg, gnu95, gnu, intelem, intel, none [gnu95]'),
        ('f2py=', None,
         'specify name (or full path) of f2py executable [{0}]'.format(
        default_f2py())),
        ]


def rebuild_static_docs(dist=None, pythondir=None):
    """Rebuild the 'static' documentation in Doc/build"""
    if pythondir:
        env = os.environ.copy()
        if 'PYTHONPATH' in env:
            env['PYTHONPATH'] = pythondir + ':' + env['PYTHONPATH']
        else:
            env['PYTHONPATH'] = pythondir
    else:
        env = None
    builddir = os.path.join(os.path.join('Doc', 'build', 'doctrees'))
    indir = os.path.join('Doc', 'source')
    outdir = os.path.join('Doc', 'build', 'html')
    cmd = 'sphinx-build -b html -d {0} {1} {2}'.format(
        builddir, indir, outdir)
    subprocess.check_call(cmd.split(), env=env)
    os.chdir('Doc')
    try:
        cmd = 'make latexpdf'
        subprocess.check_call(cmd.split(), env=env)
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

    def compile_irbempy(self):
        fcompiler = self.fcompiler
        if fcompiler == 'none':
            self.distribution.add_warning(
                'Fortran compiler specified was "none."\n'
                'IRBEM will not be available.')
            return
        # 64 bit or 32 bit?
        bit = len('%x' % sys.maxsize)*4
        irbemdir = 'irbem-lib-2012-12-12-rev425'
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
        shutil.copy(
            os.path.join(builddir, 'source', 'wrappers_{0}.inc'.format(bit)),
            os.path.join(builddir, 'source', 'wrappers.inc'.format(bit)))

        # compile irbemlib
        olddir = os.getcwd()
        os.chdir(builddir)
        F90files = ['source/onera_desp_lib.f', 'source/CoordTrans.f', 'source/AE8_AP8.f', 'source/find_foot.f']
        functions = ['make_lstar1', 'make_lstar_shell_splitting1', 'find_foot_point1',\
                     'coord_trans1','find_magequator1', 'find_mirror_point1',
                     'get_field1', 'get_ae8_ap8_flux', 'fly_in_nasa_aeap1',
                     'trace_field_line2_1', 'trace_field_line_towards_earth1']

        # call f2py
        os.system('{0} --overwrite-signature -m irbempylib -h irbempylib.pyf '
                  '{1} only: {2} :'.format(
            self.f2py, ' '.join(F90files), ' '.join(functions)))
        # intent(out) substitute list
        outlist = ['lm', 'lstar', 'blocal', 'bmin', 'xj', 'mlt', 'xout', 'bmin', 'posit', \
                   'xgeo', 'bmir', 'bl', 'bxgeo', 'flux', 'ind', 'xfoot', 'bfoot', 'bfootmag']

        inlist = ['sysaxesin', 'sysaxesout', 'iyr', 'idoy', 'secs', 'xin', 'kext', 'options', 
                  'sysaxes', 'UT', 'xIN1', 'xIN2', 'xIN3', 'stop_alt', 'hemi_flag', 'maginput']
        fln = 'irbempylib.pyf'
        if not os.path.isfile(fln):
            self.distribution.add_warning(
                'f2py failed; '
                'IRBEM will not be available.')
            os.chdir(olddir)
            return
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
            'intel': 'ifort -c -Bstatic -assume 2underscores -O2 -fPIC *.f',
            }
        compile_cmd64 = {
            'pg': 'pgf77 -c -Mnosecond_underscore -w -fastsse -fPIC *.f',
            'gnu': 'g77 -c -w -m64 -mno-align-double -O2 -fPIC -fno-second-underscore *.f',
            'gnu95': 'gfortran -m64 -c -w -O2 -fPIC -ffixed-line-length-none *.f',
            'intel': 'ifort -c -Bdynamic -O2 -fPIC *.f',
            'intelem': 'ifort -c -Bdynamic -O2 -fPIC *.f',
            }
        f2py_flags = '--fcompiler={0}'.format(fcompiler)
        if fcompiler == 'gnu':
            if bit == 32:
                f2py_flags += ' --f77flags=-fno-second-underscore,-mno-align-double'
            else:
                f2py_flags += ' --f77flags=-fno-second-underscore,-mno-align-double,-m64'
        if self.compiler:
            f2py_flags += ' --compiler={0}'.format(self.compiler)
        try:
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
        except:
            self.distribution.add_warning(
                'irbemlib compile failed. '
                'Try a different Fortran compiler? (--fcompiler)')
            os.chdir(olddir)
            return
        os.chdir('..')
        try:
            subprocess.check_call(
                '{0} -c irbempylib.pyf source/onera_desp_lib.f -Lsource -lBL2 '
                '{1}'.format(
                    self.f2py, f2py_flags),
                shell=True, env=f2py_environment(self.fcompiler))
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
        env = os.environ.copy()
        if 'PYTHONPATH' in env:
            env['PYTHONPATH'] = self.build_lib + ':' + env['PYTHONPATH']
        else:
            env['PYTHONPATH'] = self.build_lib
        cmd = 'sphinx-build -b html -d {0} {1} {2}'.format(
            builddir, indir, outdir)
        try:
            subprocess.check_call(cmd.split(), env=env)
        except:
            self.distribution.add_warning(
                "Building docs failed. Help will not be available.")

    def run(self):
        """Actually perform the build"""
        self.compile_irbempy()
        self.compile_libspacepy()
        _build.run(self)
        delete_old_files(self.build_lib)
        if self.build_docs:
            self.make_docs()
        else:
            self.copy_docs()


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
            print('Checking for dateutil...')
            import dateutil
        except:
            raise Exception("""SpacePy requires dateutil.\n
            matplotlib is the recommended way of meeing this requirement.""")
        bad = False
        try:
            print ('Checking for scipy...')
            import scipy
            print ('Checking for matplotlib...')
            import matplotlib
            assert matplotlib.compare_versions(
                matplotlib.__version__, '0.99.0')
        except:
            self.distribution.add_warning(
                'SciPy and matplotlib were not found. '
                'They are required for large parts of SpacePy.')
            bad = True
        try:
            print ('Checking for h5py...')
            import h5py
        except:
            self.distribution.add_warning(
                'h5py not found; required for parts of datamodel.')
            bad = True
        try:
            print ('Checking for networkx...')
            import networkx
        except:
            wtext = 'networkx not found; required for LANLstar.\n' + \
                    '  - see http://networkx.github.com/'
            self.distribution.add_warning(wtext)
            bad = True
        try:
            print ('Checking for ffnet...')
            import ffnet
        except:
            wtext = 'ffnet not found; required for LANLstar.\n' + \
                    '  - see http://ffnet.sourceforge.net/install.html'
            self.distribution.add_warning(wtext)
            bad = True
        if not bad:
            print('Dependencies OK.')
        _install.run(self)
        delete_old_files(self.install_lib)


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
                         if f[-4:].lower() == '.dll' and
                         (f[:11] == 'libgfortran' or
                          f[:8] == 'libgcc_s' or
                          f[:11] == 'libquadmath')]
            if len(fortnames) == 3:
                fortdir = p
                break
        if fortdir is None:
            raise RuntimeError("Can't locate fortran libraries.")
        outdir = os.path.join(self.bdist_dir, 'PLATLIB', 'spacepy', 'mingw')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for f in fortnames:
            shutil.copy(os.path.join(fortdir, f), outdir)

    def run(self):
        self.copy_fortran_libs()
        _bdist_wininst.run(self)


class sdist(_sdist):
    """Rebuild the docs before making a source distribution"""

    user_options = _sdist.user_options + compiler_options

    def initialize_options(self):
        initialize_compiler_options(self)
        _sdist.initialize_options(self)

    def finalize_options(self):
        _sdist.finalize_options(self)
        finalize_compiler_options(self)

    def run(self):
        buildcmd = self.get_finalized_command('build')
        buildcmd.run()
        rebuild_static_docs(self.distribution, buildcmd.build_lib)
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


packages = ['spacepy', 'spacepy.irbempy', 'spacepy.pycdf',
            'spacepy.plot', 'spacepy.pybats', 'spacepy.toolbox']
#If adding to package_data, also put in MANIFEST.in
package_data = ['data/*.*', 'pybats/sample_data/*', 'data/LANLstar/*']
pybats_ext = Extension('spacepy.pybats.ctrace2d',
                       sources=['spacepy/pybats/ctrace2dmodule.c'],
                       include_dirs=[numpy.get_include()])

# run setup from distutil
setup(name='spacepy',
      version='0.1.4',
      description='SpacePy: Tools for Space Science Applications',
      ext_modules=[pybats_ext, ],
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
