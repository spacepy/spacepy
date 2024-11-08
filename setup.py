#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
""" setup script to install spacepy

Authors
-------
The SpacePy Team
Los Alamos National Laboratory

Copyright 2010 - 2014 Los Alamos National Security, LLC.
"""

import sys

import setuptools
import wheel

import copy
import os, shutil, getopt, glob, re
import platform
import subprocess
if sys.platform == 'win32' and sys.version_info < (3, 8):
    # https://github.com/python/cpython/issues/84006
    from distutils import sysconfig
else:
    import sysconfig
import warnings

import distutils.ccompiler
from setuptools import setup
from distutils.command.build import build as _build
from setuptools.command.install import install as _install
from setuptools.command.build_ext import build_ext as _build_ext
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
has_editable_wheel = False
try:
    from setuptools.command.editable_wheel import editable_wheel as _editable_wheel
    has_editable_wheel = True
except ModuleNotFoundError:  # Only in very new setuptools
    pass
has_develop = False
try:
    from setuptools.command.develop import develop as _develop
    has_develop = True
except ModuleNotFoundError:  # Used in older setuptools
    pass

try:
    import setuptools.modified  # added setuptools 69.0.0
except ImportError:
    import setuptools.dep_util  # removed setuptools 70.0.0
    setuptools.modified = setuptools.dep_util
import setuptools.extension

import distutils.sysconfig


# building official release, fail fast instead of accepting "partly works"
release_build = bool(os.environ.get("SPACEPY_RELEASE", False))

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


def initialize_compiler_options(cmd):
    """Initialize the compiler options for a command"""
    cmd.compiler = None


def finalize_compiler_options(cmd):
    """Finalize compiler options for a command

    If compiler options (compiler) have not been
    specified for a command, check if they were specified for other
    commands on the command line--if so, grab from there. If not,
    set reasonable defaults.

    cmd: the command to finalize the options for
    """
    dist = cmd.distribution
    defaults = {'compiler': None,
                }
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
    if cmd.compiler == None and sys.platform == 'win32':
        cmd.compiler = 'mingw32'


class build(_build):
    """Support compiler options on build"""

    def initialize_options(self):
        _build.initialize_options(self)
        initialize_compiler_options(self)

    def finalize_options(self):
        _build.finalize_options(self)
        finalize_compiler_options(self)


class build_ext(_build_ext):
    """Extends base distutils build_ext to make libspacepy, irbem"""

    def initialize_options(self):
        _build_ext.initialize_options(self)
        initialize_compiler_options(self)

    def finalize_options(self):
        _build_ext.finalize_options(self)
        finalize_compiler_options(self)

    def compile_irbempy(self):
        """Compile the irbem library

        Returns path to compiled shared library if successful.
        """
        # 64 bit or 32 bit?
        bit = len('%x' % sys.maxsize)*4
        irbemdir = 'irbem-lib-20220829-dfb9d26'
        srcdir = os.path.join('spacepy', 'irbempy', irbemdir, 'source')
        outdir = os.path.join(os.path.abspath(self.build_lib),
                              'spacepy', 'irbempy')

        # Check for existing shared library and if up to date
        ccomp = distutils.ccompiler.new_compiler(compiler=self.compiler)
        if hasattr(distutils.ccompiler, 'customize_compiler'):
            distutils.ccompiler.customize_compiler(ccomp)
        else:
            distutils.sysconfig.customize_compiler(ccomp)
        irbemname = ccomp.library_filename('irbem', lib_type='shared')
        libirbem = os.path.join(outdir, irbemname)
        if os.path.exists(libirbem):
            sources = glob.glob(os.path.join(srcdir, '*.f')) + \
                      glob.glob(os.path.join(srcdir, '*.inc'))
            if not setuptools.modified.newer_group(sources, libirbem):
                return irbemname

        if not sys.platform in ('darwin', 'linux2', 'linux', 'win32'):
            warnings.warn(
                '%s not supported at this time. ' % sys.platform +
                'IRBEM will not be available')
            return None

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
        print('Building irbem library...')
        # compile (platform dependent)
        os.chdir('source')
        compflags = ['-w', '-O2', '-fPIC', '-ffixed-line-length-none',
                     '-std=legacy']
        if not sys.platform.startswith('win') \
           and not platform.uname()[4].startswith(('arm', 'aarch64')):
            # Raspberry Pi doesn't have or need this switch
            compflags = ['-m{0}'.format(bit)] + compflags
        comp_candidates = ['gfortran']
        if sys.platform == 'win32':
            comp_candidates.insert(0, 'gfortran.exe')
        if 'FC' in os.environ:
            comp_candidates.insert(0, os.environ['FC'])
        for fc in comp_candidates:
            retval = subprocess.call([fc, '-c'] + compflags
                                     + list(glob.glob('*.f')))
            if retval == 0:
                break
            else:
                warnings.warn('Compiler {0} failed, trying another'.format(fc))
        else:
            if release_build:
                raise RuntimeError('irbemlib compile failed.')
            warnings.warn('irbemlib compile failed.')
            os.chdir(olddir)
            return None
        # fc is known-good compiler, link shared object
        ldflags = os.environ.get('LDFLAGS', '').split()
        if '-shared' not in ldflags:
            ldflags.insert(0, '-shared')
        if '-fPIC' in ldflags:
            del ldflags[ldflags.index('-fPIC')]
        if sys.platform == 'darwin':
            if not any([flag.startswith(f'-mmacosx-version-min') for flag in ldflags]):
                isarm = platform.uname()[4].startswith(('arm', 'aarch64'))
                min_os_ver = 11.0 if isarm else 10.9
                ldflags.insert(1, f"-mmacosx-version-min={min_os_ver}")
            if 'SDKROOT' in os.environ and '-isysroot' not in ldflags:
                ldflags.insert(1, '-isysroot')
                ldflags.insert(2, os.environ['SDKROOT'])
        link_irbemlib = [fc] + ldflags \
            + list(glob.glob('*.o')) \
            + ['-o', irbemname, '-fPIC']
        subprocess.check_call(link_irbemlib)
        if not os.path.exists(irbemname):
            if release_build:
                raise RuntimeError('irbem shared library failed.')
            warnings.warn(
                'irbemlib shared library failed.')
            os.chdir(olddir)
            return None
        shutil.move(irbemname, libirbem)

        if sys.platform == 'darwin':
            # Look for the library location that is shipped with the wheel
            cmd = ['install_name_tool', '-add_rpath', '@loader_path/../libs/',
                   libirbem]
            try:
                subprocess.call(cmd)
            except FileNotFoundError:
                if release_build:
                    raise
        os.chdir(olddir)
        return libirbem

    def compile_libspacepy(self):
        """Compile the C library, libspacepy

        Returns path to the library if successful
        """
        srcdir = os.path.join('spacepy', 'libspacepy')
        outdir = os.path.join(self.build_lib, 'spacepy')
        try:
            comp = distutils.ccompiler.new_compiler(compiler=self.compiler)
            if sys.platform == 'win32':
                # Cut out MSVC runtime https://bugs.python.org/issue16472
                comp.dll_libraries = []
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
                        if setuptools.modified.newer_group([s] + headers, o)]
            if outdated:
                comp.compile(outdated, output_dir=self.build_temp)
            libpath = os.path.join(
                outdir, comp.library_filename('spacepy', lib_type='shared'))
            if setuptools.modified.newer_group(objects, libpath):
                comp.link_shared_lib(objects, 'spacepy', libraries=['m'],
                                     output_dir=outdir)
            return libpath
        except:
            if release_build:
                raise
            warnings.warn(
                'libspacepy compile failed; some operations may be slow.')
            print('libspacepy compile failed:')
            (t, v, tb) = sys.exc_info()
            print(v)

    def run(self):
        """Actually perform the extension build"""
        libspacepy = self.compile_libspacepy()
        libirbem = self.compile_irbempy()
        self._outputs = [l for l in (libspacepy, libirbem)
                         if l is not None]
        if sys.platform == 'win32':
            #Copy mingw32 DLLs. This keeps them around if ming is uninstalled,
            #but more important puts them where bdist_wheel
            #will include them in binary installers
            libs = copy_win_libs(os.path.join(self.build_lib, 'spacepy'))
            self._outputs.extend(libs)
        if sys.platform == 'darwin' and libirbem:
            # Copy gfortran dyanamic libraries
            # Puts them where bdist_wheel will include them in binary installers
            libs = copy_mac_libs(os.path.join(self.build_lib, 'spacepy'))
            self._outputs.extend(libs)
        if sys.platform == 'linux':
            # Copy CDF library where bdist_wheel will include it in binary
            libs = copy_linux_libs(os.path.join(self.build_lib, 'spacepy'))
            self._outputs.extend(libs)
        if not (getattr(self, 'editable_mode', False)
                or getattr(self, 'inplace', False)):
            return
        # Copy compiled outputs into the source
        build_py = self.distribution.get_command_obj('build_py')
        package_dir = build_py.get_package_dir('spacepy')
        if libspacepy is not None and os.path.exists(libspacepy):
            shutil.copy2(libspacepy, package_dir)
        if libirbem is not None and os.path.exists(libirbem):
            shutil.copy2(libirbem, os.path.join(package_dir, 'irbempy'))

    def get_outputs(self):
        return self._outputs


class install(_install):
    """Support compiler options on install"""

    user_options = _install.user_options

    def initialize_options(self):
        initialize_compiler_options(self)
        _install.initialize_options(self)

    def finalize_options(self):
        _install.finalize_options(self)
        finalize_compiler_options(self)


def copy_win_libs(outdir):
    """Copy pre-built (binary) Windows libraries into a build

    This includes the mingw runtime and CDF libraries

    :param str outdir: Final target directory of the DLLs in the build.
    :returns list: List of copied libraries
    """
    outputs = []
    libdir = None
    libnames = None
    libneeded = ('libgfortran', 'libgcc_s', 'libquadmath',)
    liboptional = ('libwinpthread',)
    for p in os.environ['PATH'].split(';'):
        if not os.path.isdir(p):
            continue
        libnames = [
            f for f in os.listdir(p) if f.lower().endswith('.dll')
            and f.startswith(libneeded + liboptional)]
        if len([f for f in libnames if f.startswith(libneeded)])\
           == len(libneeded):
            libdir = p
            break
    if libdir is None:
        if release_build:
            raise RuntimeError("Can't locate runtime libraries.")
        return []
    mingdir = os.path.join(outdir, 'libs')
    if not os.path.exists(mingdir):
        os.makedirs(mingdir, exist_ok=True)
    for f in libnames:
        shutil.copy(os.path.join(libdir, f), mingdir)
    outputs.extend((os.path.join(mingdir, f) for f in libnames))
    # Copy CDF library from root spacepy directory
    if isinstance(__file__, str):
        cdfdll = os.path.join(os.path.dirname(__file__), 'dllcdf.dll')
        if os.path.isfile(cdfdll):
            shutil.copy(cdfdll, outdir)
            outputs.append(os.path.join(outdir, 'dllcdf.dll'))
        elif release_build:
            raise RuntimeError(f'{cdfdll} is not a file')
    elif release_build:
        raise RuntimeError('Unable to find SpacePy source to copy CDF library.')
    return outputs


def copy_mac_libs(outdir):
    """Copy Mac fortran and CDF libraries into a build

    :param str outdir: Final target directory of the libraries in the build.
    :returns list: List of copied libraries
    """
    outputs = []
    libneeded = ['libgfortran', 'libgcc_s.1.1', 'libgcc_s.1', 'libquadmath']
    outlibdir = os.path.join(outdir, 'libs')
    if not os.path.exists(outlibdir):
        os.makedirs(outlibdir, exist_ok=True)
    for f in libneeded:
        p = subprocess.Popen(['gfortran', f'--print-file-name={f}.dylib'],
                         stdout=subprocess.PIPE)
        out, err = p.communicate()
        if p.returncode:
            if release_build:
                raise RuntimeError(f"Failure searching for {f}.")
            return []
        libfile = out.rstrip().decode('ascii')
        if not os.path.isfile(libfile):
            if f.startswith('libgcc_s'):
                continue  # Only need one of them
            if release_build:
                raise RuntimeError(f"{f} not found.")
            return []
        libfile = os.path.realpath(libfile)
        shutil.copy(libfile, outlibdir)
        outputs.append(os.path.join(outlibdir, os.path.basename(libfile)))
        if f != 'libgcc_s.1':
            continue
        # In some cases libgcc_s.1 has a hardcoded ref to 1.1 that needs fixing
        libgcc_s = outputs[-1]
        p = subprocess.Popen(['otool', '-L', libgcc_s],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if p.returncode:
            if release_build:
                raise RuntimeError("Failed to read library paths from libgcc_s")
            return []
        # libgcc_s.1 refers to .1.1 with hardcoded path, fix that
        hardcoded = [l.strip().split()[0]
                     for l in out.decode('ascii').split('\n')
                     if l.startswith('\t/') and 'libgcc_s.1.1' in l]
        if not hardcoded:  # No reference, no worries
            continue
        assert len(hardcoded) == 1
        cmd = ['install_name_tool', '-change', hardcoded[0],
               os.path.join('@loader_path', os.path.basename(hardcoded[0])),
               libgcc_s]
        try:
            subprocess.call(cmd)
        except FileNotFoundError:
            if release_build:
                raise
    if isinstance(__file__, str):
        cdflib = os.path.join(os.path.dirname(__file__), 'libcdf.dylib')
        if os.path.isfile(cdflib):
            shutil.copy(cdflib, outdir)
            outputs.append(os.path.join(outdir, 'libcdf.dylib'))
        elif release_build:
            raise RuntimeError(f'{cdflib} is not a file')
    elif release_build:
        raise RuntimeError('Unable to find SpacePy source to copy CDF library.')
    if release_build and not [f for f in outputs if 'libgcc' in f]:
        raise RuntimeError('No libgcc found.')
    return outputs


def copy_linux_libs(outdir):
    """Copy Linux CDF libraries into a build

    :param str outdir: Final target directory of the libraries in the build.
    :returns list: List of copied libraries
    """
    outputs = []
    if isinstance(__file__, str):
        cdflib = os.path.join(os.path.dirname(__file__), 'libcdf.so')
        if os.path.isfile(cdflib):
            shutil.copy(cdflib, outdir)
            outputs.append(os.path.join(outdir, 'libcdf.so'))
        elif release_build:
            raise RuntimeError(f'{cdflib} is not a file')
    elif release_build:
        raise RuntimeError('Unable to find SpacePy source to copy CDF library.')
    return outputs


class bdist_wheel(_bdist_wheel):
    """Handle compiler options for wheel build"""

    def initialize_options(self):
        initialize_compiler_options(self)
        _bdist_wheel.initialize_options(self)

    def finalize_options(self):
        _bdist_wheel.finalize_options(self)
        finalize_compiler_options(self)


if has_editable_wheel:
    class editable_wheel(_editable_wheel):
        """Handle compiler options for editable wheel build"""

        def initialize_options(self):
            initialize_compiler_options(self)
            _editable_wheel.initialize_options(self)

        def finalize_options(self):
            _editable_wheel.finalize_options(self)
            finalize_compiler_options(self)


if has_develop:
    class develop(_develop):
        """Make sure old-style editable install has compiler options"""

        def initialize_options(self):
            initialize_compiler_options(self)
            _develop.initialize_options(self)

        def finalize_options(self):
            _develop.finalize_options(self)
            finalize_compiler_options(self)


packages = ['spacepy', 'spacepy.irbempy', 'spacepy.pycdf',
            'spacepy.plot', 'spacepy.pybats', 'spacepy.toolbox',
            'spacepy.ctrans',
            'spacepy.data', 'spacepy.data.LANLstar',
            'spacepy.data.TS07D.TAIL_PAR',
            ]
#If adding to package_data, also put in MANIFEST.in
package_data = {
    'spacepy.data': ['*.*'],
    'spacepy.data.LANLstar': ['*'],
    'spacepy.data.TS07D.TAIL_PAR': ['*'],
}
# Built with custom code that handles the source files
ext_modules = [setuptools.extension.Extension('spacepy.irbempy.irbempylib', [])]

# Duplicated between here and pyproject.toml because pyproject.toml support
# requires setuptools 61.0.0, and doesn't support extensions
setup_kwargs = {
    'name': 'spacepy',
    'version': '0.8.0a0',
    'description': 'SpacePy: Tools for Space Science Applications',
    'long_description': 'SpacePy: Tools for Space Science Applications',
    'author': 'SpacePy team',
    'author_email': 'spacepy@lanl.gov',
    'maintainer': 'Steve Morley, Dan Welling, Brian Larsen, Jon Niehof',
    'maintainer_email': 'spacepy@lanl.gov',
    'url': 'https://github.com/spacepy/spacepy',
    'packages': packages,
    'package_data': package_data,
    'ext_modules': ext_modules,
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Python Software Foundation License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: C',
        'Programming Language :: Fortran',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    'keywords': ['magnetosphere', 'plasma', 'physics', 'space', 'solar.wind', 'space.weather', 'magnetohydrodynamics'],
    'license':  'PSF',
    'platforms':  ['Windows', 'Linux', 'MacOS X', 'Unix'],
    'install_requires': [
        'numpy>=1.15.1',
        'scipy>=1.0',
        'matplotlib>=3.1',
        'h5py>=2.10',
        'python_dateutil>=2.5',
        # AstroPy is only required to convert to/from AstroPy, so either
        # user has it or don't care.
        #'astropy>=2.0',
        # Similar for pandas
        #'pandas>=0.18',
    ],
    'python_requires': '>=3.7',
    'cmdclass': {'build': build,
                 'build_ext': build_ext,
                 'install': install,
                 'bdist_wheel': bdist_wheel,
          },
    'zip_safe': False,
    'options': {'bdist_wheel': {'py_limited_api': 'cp37'}},
}

if has_editable_wheel:
    setup_kwargs['cmdclass']['editable_wheel'] = editable_wheel
if has_develop:
    setup_kwargs['cmdclass']['develop'] = develop

# run setup from distutil
with warnings.catch_warnings(record=True) as warnlist:
    setup(**setup_kwargs)
    if len(warnlist):
        print('\nsetup produced the following warnings. '
              'Some functionality may be missing.\n')
        for w in warnlist:
            print(w.message)
