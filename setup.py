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

#pip force-imports setuptools, on INSTALL, so then need to use its versions
#but on reading the egg info, it DOESN'T force-import, assumes you are using
import sys
if any([a in sys.argv for a in ('pip-egg-info', 'bdist_wheel')]):
    import setuptools
if 'bdist_wheel' in sys.argv:
    import wheel
use_setuptools = "setuptools" in globals()
#Calling egg info, so a lot of stuff doesn't have to work
egginfo_only = ('pip-egg-info' in sys.argv)

import copy
import os, shutil, getopt, glob, re
import subprocess
import warnings

import distutils.ccompiler
#Save these because numpy will smash them, but need numpy for Fortran stuff
#numpy REPLACES new_compiler function, and EDITS compiler class dict
real_compiler_class = copy.deepcopy(distutils.ccompiler.compiler_class)
real_distutils_ccompiler_new_compiler = distutils.ccompiler.new_compiler
try:
    from numpy.distutils.core import setup
    from numpy.distutils.command.build import build as _build
    #The numpy versions automatically use setuptools if necessary
    #We need the numpy setup to handle fortran compiler options
    from numpy.distutils.command.install import install as _install
    from numpy.distutils.command.sdist import sdist as _sdist
except: #numpy not installed, hopefully just getting egg info
    if not egginfo_only:
        raise
    from distutils.core import setup
    from distutils.command.build import build as _build
    if use_setuptools:
        from setuptools.command.install import install as _install
        from setuptools.command.sdist import sdist as _sdist
    else:
        from distutils.command.install import install as _install
        from distutils.command.sdist import sdist as _sdist

if use_setuptools:
    from setuptools.command.bdist_wininst import bdist_wininst as _bdist_wininst
else:
    from distutils.command.bdist_wininst import bdist_wininst as _bdist_wininst
if 'bdist_wheel' in sys.argv:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
import distutils.ccompiler
import distutils.dep_util
import distutils.sysconfig
from distutils.errors import DistutilsOptionError
try:
    import importlib.machinery #Py3.3 and later
except ImportError:
    import imp #pre-3.3
else:
    #importlib.machinery exists in 3.2, but doesn't have this
    if hasattr(importlib.machinery, 'ExtensionFileLoader'):
        imp = None
    else:
        import imp #fall back to old-style
try:
    import numpy
    import numpy.distutils.command.config_compiler
except:
    if egginfo_only:
        pass
    else:
        raise


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
               'pycdf/_pycdf.*', 'toolbox/toolbox.py*', ]


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


def f2py_options(fcompiler, dist=None):
    """Get an OS environment for f2py, and find name of Fortan compiler

    The OS environment puts in the shared options if LDFLAGS is set
    """
    env = None
    import numpy.distutils.fcompiler
    numpy.distutils.fcompiler.load_all_fcompiler_classes()
    if not fcompiler in numpy.distutils.fcompiler.fcompiler_class:
        return (None, None) #and hope for the best
    fcomp = numpy.distutils.fcompiler.fcompiler_class[fcompiler][1]
    #Various compilers specify executables for ranlib/ar, but the base
    #class command_vars explicitly ignores them. So monkeypatch.
    for k in ('archiver', 'ranlib'):
        if k in fcomp.executables and k in fcomp.command_vars._conf_keys:
            oldval = fcomp.command_vars._conf_keys[k]
            if oldval[0] is None:
                fcomp.command_vars._conf_keys[k] = ('exe.{0}'.format(k),) +  \
                                                   oldval[1:]
    fcomp = fcomp()
    fcomp.customize(dist)
    if 'LDFLAGS' in os.environ:
        env = os.environ.copy()
        currflags = os.environ['LDFLAGS'].split()
        fcompflags = fcomp.get_flags_linker_so()
        it = iter(range(len(fcompflags)))
        for i in it:
            if i == len(fcompflags) - 1 or fcompflags[i + 1].startswith('-'):
                #a simple flag
                if not fcompflags[i] in currflags:
                    currflags.append(fcompflags[i])
                continue
            #Flag that takes an option, consume the option (and maybe add)
            next(it)
            idx = 0
            while fcompflags[i] in currflags[idx:]:
                idx = currflags.index(fcompflags[i], idx)
                if idx < len(currflags) + 1 and \
                   currflags[idx + 1] == fcompflags[i + 1]:
                    break
            else:
                #Was NOT found, so add it
                currflags.append(fcompflags[i])
                currflags.append(fcompflags[i + 1])
        env['LDFLAGS'] = ' '.join(currflags)
    return (env, fcomp.executables)


def initialize_compiler_options(cmd):
    """Initialize the compiler options for a command"""
    cmd.fcompiler = None
    cmd.f2py = None
    cmd.compiler = None
    cmd.f77exec = None
    cmd.f90exec = None


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
                'compiler': None,
                'f77exec': None,
                'f90exec': None,}
    optdict = dist.get_option_dict(cmd.get_command_name())
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
        #Also explicitly carry over command line option, needed for config_fc
        if not option in optdict:
            for c in dist.commands:
                otheroptdict = dist.get_option_dict(c)
                if option in otheroptdict:
                    optdict[option] = otheroptdict[option]
                    break
    #Special-case defaults, checks
    if not cmd.fcompiler in ('pg', 'gnu', 'gnu95', 'intelem', 'intel', 'none', 'None'):
        raise DistutilsOptionError(
            '--fcompiler={0} unknown'.format(cmd.fcompiler) +
            ', options: pg, gnu, gnu95, intelem, intel, None')
    if len('%x' % sys.maxsize)*4 == 32 and cmd.fcompiler == 'intelem':
        raise DistutilsOptionError(
            '--fcompiler=intelem requires a 64-bit architecture')
    if cmd.compiler == None and sys.platform == 'win32':
        cmd.compiler = 'mingw32'
    #Add interpreter to f2py if it needs it (usually on Windows)
    #If it's a list, it's already been patched up
    if sys.platform == 'win32' and isinstance(cmd.f2py, str) \
       and not is_win_exec(cmd.f2py):
        if egginfo_only: #Punt, we're not going to call it
            cmd.f2py = [cmd.f2py]
            return
        f2py = cmd.f2py
        if not os.path.isfile(f2py): #Not a file, and didn't exec
            f2pydir = next((d for d in os.environ['PATH'].split(os.pathsep)
                            if os.path.isfile(os.path.join(d, f2py))), None)
            if f2pydir: #Found the file, try it
                f2py = os.path.join(f2pydir, f2py)
                if not is_win_exec(f2py): #Try the interpreter
                    if is_win_exec(sys.executable, f2py):
                        cmd.f2py = [sys.executable, f2py]
                    else: #Nothing to be done
                        raise RuntimeError(
                            'f2py {0} found but not executable'.format(
                                cmd.f2py))
                else: #Found and executable (unlikely, since would have worked)
                    cmd.f2py = [f2py]
            else: #Couldn't find the file, just try the interpreter
                if is_win_exec(sys.executable, f2py):
                    cmd.f2py = [sys.executable, f2py]
                else: #Nothing to be done
                    raise RuntimeError(
                        'f2py {0} not found and not executable'.format(
                            cmd.f2py))
        else: #Not executable, but IS a file
            if is_win_exec(sys.executable, f2py):
                cmd.f2py = [sys.executable, f2py]
            else: #Nothing to be done
                raise RuntimeError(
                    'f2py {0} exists but not executable'.format(
                        cmd.f2py))
    else:
        cmd.f2py = [cmd.f2py]


def is_win_exec(*args):
    """Test if a file spec is an executable

    This is really only useful on Windows

    :param list args: arguments to call
    :returns: True if the arguments can be called with subprocess
    :rtype: bool
    """
    try:
        p = subprocess.Popen(args, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    except WindowsError:
        return False
    else:
        output, errors = p.communicate()
        return not(p.returncode)


compiler_options = [
        ('fcompiler=', None,
         'specify the fortran compiler to use: pg, gnu95, gnu, intelem, intel, none [gnu95]'),
        ('f2py=', None,
         'specify name (or full path) of f2py executable [{0}]'.format(
        default_f2py())),
        ('f77exec=', None,
         'specify the path to the F77 compiler'),
        ('f90exec=', None,
         'specify the path to the F90 compiler'),
        ]


def rebuild_static_docs():
    """Rebuild the 'static' documentation in Doc/build"""
    builddir = os.path.join(os.path.join('Doc', 'build', 'doctrees'))
    indir = os.path.join('Doc', 'source')
    outdir = os.path.join('Doc', 'build', 'html')
    cmd = '{0} -b html -d {1} {2} {3}'.format(
        os.environ['SPHINXBUILD'] if 'SPHINXBUILD' in os.environ
        else 'sphinx-build',
        builddir, indir, outdir)
    subprocess.check_call(cmd.split())
    os.chdir('Doc')
    try:
        cmd = '{0}{1} latexpdf'.format(
            os.environ['MAKE'] if 'MAKE' in os.environ else 'make',
            ('SPHINXBUILD=' + os.environ['SPHINXBUILD'])
            if 'SPHINXBUILD' in os.environ else '')
        subprocess.check_call(cmd.split())
    except:
        warnings.warn('PDF documentation rebuild failed:')
        (t, v, tb) = sys.exc_info()
        print(v)
    finally:
        os.chdir('..')


#Possible names of the irbem output library. Unfortunately this seems
#to depend on Python version, f2py version, and phase of the moon
def get_irbem_libfiles():
    cvars = distutils.sysconfig.get_config_vars()
    libfiles = ['irbempylib' + cvars[ext] for ext in ('SO', 'EXT_SUFFIX')
                if ext in cvars]
    if len(libfiles) < 2: #did we get just the ABI-versioned one?
        abi = distutils.sysconfig.get_config_var('SOABI')
        if abi and libfiles[0].startswith('irbempylib.' + abi):
            libfiles.append('irbempylib' +
                            libfiles[0][(len('irbempylib.') + len(abi)):])
    if len(libfiles) == 2 and libfiles[0] == libfiles[1]:
        del libfiles[0]
    return libfiles


class build(_build):
    """Extends base distutils build to make pybats, libspacepy, irbem"""

    if not egginfo_only:
        sub_commands = [('config_fc', lambda *args:True)] + _build.sub_commands

    user_options = _build.user_options + compiler_options + [
        ('build-docs', None,
         'Build documentation with Sphinx (default: copy pre-built) [False]'),
        ]

    def initialize_options(self):
        self.build_docs = None
        _build.initialize_options(self)
        initialize_compiler_options(self)

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
        if fcompiler in ['none', 'None']:
            warnings.warn(
                'Fortran compiler specified was "none."\n'
                'IRBEM will not be available.')
            return
        # 64 bit or 32 bit?
        bit = len('%x' % sys.maxsize)*4
        irbemdir = 'irbem-lib-2019-04-04-rev620'
        srcdir = os.path.join('spacepy', 'irbempy', irbemdir, 'source')
        outdir = os.path.join(os.path.abspath(self.build_lib),
                              'spacepy', 'irbempy')
        #Possible names of the output library. Unfortunately this seems
        #to depend on Python version, f2py version, and phase of the moon
        libfiles = get_irbem_libfiles()
        #Delete any irbem extension modules from other versions
        for f in glob.glob(os.path.join(outdir, 'irbempylib*')):
            if not os.path.basename(f) in libfiles:
                os.remove(f)
        #Does a matching one exist?
        existing_libfiles = [f for f in libfiles
                             if os.path.exists(os.path.join(outdir, f))]
        #Can we import it?
        importable = []
        for f in existing_libfiles:
            fspec = os.path.join(outdir, f)
            if imp: #old-style imports
                suffixes = imp.get_suffixes()
                desc = next(
                    (s for s in imp.get_suffixes() if f.endswith(s[0])), None)
                if not desc: #apparently not loadable
                    os.remove(fspec)
                    continue
                fp = open(fspec, 'rb')
                try:
                    imp.load_module('irbempylib', fp, fspec, desc)
                except ImportError:
                    fp.close()
                    os.remove(fspec)
                else:
                    fp.close()
                    importable.append(f)
            else: #Py3.3 and later imports, not tested
                loader = importlib.machinery.ExtensionFileLoader(
                    'irbempylib', fspec)
                try:
                    loader.load_module('irbempylib')
                except ImportError:
                    os.remove(fspec)
                else:
                    importable.append(f)
        existing_libfiles = importable
        #if MORE THAN ONE matching output library file, delete all;
        #no way of knowing which is the correct one or if it's up to date
        if len(existing_libfiles) > 1:
            for f in existing_libfiles:
                os.remove(os.path.join(outdir, f))
            existing_libfiles = []
        #If there's still one left, check up to date
        if existing_libfiles:
            sources = glob.glob(os.path.join(srcdir, '*.f')) + \
                      glob.glob(os.path.join(srcdir, '*.inc'))
            if not distutils.dep_util.newer_group(
                sources, os.path.join(outdir, existing_libfiles[0])):
                return

        if not sys.platform in ('darwin', 'linux2', 'linux', 'win32'):
            warnings.warn(
                '%s not supported at this time. ' % sys.platform +
                'IRBEM will not be available')
            return
        if fcompiler == 'pg' and sys.platform == 'darwin':
            warnings.warn(
                'Portland Group compiler "pg" not supported on Mac OS\n'
                'IRBEM will not be available.')
            return
        if fcompiler != 'gnu95' and sys.platform == 'win32':
            warnings.warn(
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

        f2py_env, fcompexec = f2py_options(fcompiler, self.distribution)

        # compile irbemlib
        olddir = os.getcwd()
        os.chdir(builddir)
        F90files = ['source/onera_desp_lib.f', 'source/CoordTrans.f', 'source/AE8_AP8.f', 'source/find_foot.f',\
                    'source/LAndI2Lstar.f', 'source/drift_bounce_orbit.f']
        functions = ['make_lstar1', 'make_lstar_shell_splitting1', 'find_foot_point1',\
                     'coord_trans1','find_magequator1', 'find_mirror_point1',
                     'get_field1', 'get_ae8_ap8_flux', 'fly_in_nasa_aeap1',
                     'trace_field_line2_1', 'trace_field_line_towards_earth1', 'trace_drift_bounce_orbit',
                     'landi2lstar1', 'landi2lstar_shell_splitting1']

        # call f2py
        cmd = self.f2py + ['--overwrite-signature', '-m', 'irbempylib', '-h',
               'irbempylib.pyf'] + F90files + ['only:'] + functions + [':']
        subprocess.check_call(cmd)
        # intent(out) substitute list
        outlist = ['lm', 'lstar', 'blocal', 'bmin', 'xj', 'mlt', 'xout', 'bmin', 'posit', \
                   'xgeo', 'bmir', 'bl', 'bxgeo', 'flux', 'ind', 'xfoot', 'bfoot', 'bfootmag',\
                   'leI0', 'Bposit', 'Nposit', 'hmin', 'hmin_lon']

        inlist = ['sysaxesin', 'sysaxesout', 'iyr', 'idoy', 'secs', 'xin', 'kext', 'options', 
                  'sysaxes', 'UT', 'xIN1', 'xIN2', 'xIN3', 'stop_alt', 'hemi_flag', 'maginput',\
                  't_resol', 'r_resol', 'lati', 'longi', 'alti', 'R0','xx0']
        fln = 'irbempylib.pyf'
        if not os.path.isfile(fln):
            warnings.warn(
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

        print('Building irbem library...')
        # compile (platform dependent)
        os.chdir('source')
        comppath = {
            'pg': 'pgf77',
            'gnu': 'g77',
            'gnu95': 'gfortran',
            'intel': 'ifort',
            'intelem': 'ifort',
            }[fcompiler]
        compflags = {
            'pg': '-Mnosecond_underscore -w -fastsse -fPIC',
            'gnu': '-w -O2 -fPIC -fno-second-underscore',
            'gnu95': '-w -O2 -fPIC -ffixed-line-length-none -std=legacy',
            'intel': '-Bstatic -assume 2underscores -O2 -fPIC',
            'intelem': '-Bdynamic -O2 -fPIC',
            }[fcompiler]
        if fcompiler == 'gnu':
            if bit == 64:
                compflags = '-m64 ' + compflags
        if fcompiler == 'gnu95' and not os.uname()[4].startswith('arm'):
            # Raspberry Pi doesn't have this switch and assumes 32-bit
            compflags = '-m{0} '.format(bit) + compflags
        if fcompiler.startswith('intel'):
            if bit == 32:
                compflags = '-Bstatic -assume 2underscores ' + compflags
            else:
                compflags = '-Bdynamic ' + compflags
        comp_candidates = [comppath]
        if fcompexec is not None and 'compiler_f77' in fcompexec:
            comp_candidates.insert(0, fcompexec['compiler_f77'][0])
        for fc in comp_candidates:
            retval = subprocess.call(fc + ' -c ' + compflags + ' *.f',
                                     shell=True)
            if retval == 0:
                break
            else:
                warnings.warn('Compiler {0} failed, trying another'.format(fc))
        else:
            warnings.warn('irbemlib compile failed. '
                          'Try a different Fortran compiler? (--fcompiler)')
            os.chdir(olddir)
            return
        retval = -1
        if 'archiver' in fcompexec:
            archiver = ' '.join(fcompexec['archiver']) + ' '
            ranlib = None
            if 'ranlib' in fcompexec:
                ranlib = ' '.join(fcompexec['ranlib']) + ' '
            retval = subprocess.check_call(archiver + 'libBL2.a *.o', shell=True)
            if (retval == 0) and ranlib:
                retval = subprocess.call(ranlib + 'libBL2.a', shell=True)
            if retval != 0:
                warnings.warn(
                    'irbemlib linking failed, trying with default linker.')
        if retval != 0: #Try again with defaults
            archiver = {
                'darwin': 'libtool -static -o ',
                'linux': 'ar -r ',
                'linux2': 'ar -r ',
                'win32': 'ar - r',
                }[sys.platform]
            ranlib = {
                'darwin': None,
                'linux': 'ranlib ',
                'linux2': 'ranlib ',
                'win32': 'ranlib ',
                }[sys.platform]
            try:
                subprocess.check_call(archiver + 'libBL2.a *.o', shell=True)
                if ranlib:
                    subprocess.check_call(ranlib + 'libBL2.a', shell=True)
            except:
                warnings.warn(
                    'irbemlib linking failed. '
                    'Try a different Fortran compiler? (--fcompiler)')
                os.chdir(olddir)
                return
        os.chdir('..')

        f2py_flags = ['--fcompiler={0}'.format(fcompiler)]
        if fcompiler == 'gnu':
            f2py_flags.append('--f77flags=-fno-second-underscore,-mno-align-double')
            if bit == 64:
                f2py_flags[-1] += ',-m64'
        if fcompiler == 'gnu95':
            f2py_flags.extend(['--f77flags=-std=legacy',
                               '--f90flags=-std=legacy'])
        if self.compiler:
            f2py_flags.append('--compiler={0}'.format(self.compiler))
        if self.f77exec:
            f2py_flags.append('--f77exec={0}'.format(self.f77exec))
        if self.f90exec:
            f2py_flags.append('--f90exec={0}'.format(self.f90exec))
        try:
            subprocess.check_call(
                self.f2py + ['-c', 'irbempylib.pyf', 'source/onera_desp_lib.f',
                 '-Lsource', '-lBL2'] + f2py_flags, env=f2py_env)
        except:
            warnings.warn(
                'irbemlib module failed. '
                'Try a different Fortran compiler? (--fcompiler)')

        #All matching outputs
        created_libfiles = [f for f in libfiles if os.path.exists(f)]
        if len(created_libfiles) == 0: #no matches
            warnings.warn(
                'irbemlib build produced no recognizable module. '
                'Try a different Fortran compiler? (--fcompiler)')
        elif len(created_libfiles) == 1: #only one, no ambiguity
            shutil.move(created_libfiles[0],
                        os.path.join(outdir, created_libfiles[0]))
        elif len(created_libfiles) == 2 and \
                len(existing_libfiles) == 1: #two, so one is old and one new
            for f in created_libfiles:
                if f == existing_libfiles[0]: #delete the old one
                    os.remove(f)
                else: #and move the new one to its place in build
                    shutil.move(f,
                                os.path.join(outdir, f))
        else:
             warnings.warn(
                'irbem build failed: multiple build outputs ({0}).'.format(
                     ', '.join(created_libfiles)))
        os.chdir(olddir)
        return

    def compile_libspacepy(self):
        """Compile the C library, libspacepy. JTN 20110224"""
        srcdir = os.path.join('spacepy', 'libspacepy')
        outdir = os.path.join(self.build_lib, 'spacepy')
        try:
            if sys.platform == 'win32':
                #numpy whacks our C compiler options. We need it for Fortran,
                #but since we're using C to build a standard shared object,
                #not a numpy extension module, need vanilla C back.
                numpy_compiler_class = distutils.ccompiler.compiler_class
                distutils.ccompiler.compiler_class = real_compiler_class
                comp = real_distutils_ccompiler_new_compiler(
                    compiler=self.compiler)
                distutils.ccompiler.compiler_class = numpy_compiler_class
                #Cut out MSVC runtime https://bugs.python.org/issue16472
                #For some reason they're named differently on py2 and py3
                comp.dll_libraries = [
                    l for l in comp.dll_libraries if not l.startswith((
                        'msvcr', 'vcruntime'))]
            else:
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
            warnings.warn(
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
                warnings.warn(
                "Numpydoc and sphinx required to build documentation.\n"
                "Help will not be available; try without --build-docs.")
                return
            else:
                warnings.warn(
                "Numpydoc and sphinx required to build documentation.\n"
                "Help will not be available.")
                return
        builddir = os.path.join(os.path.join(self.build_temp, 'doctrees'))
        indir = os.path.join('Doc', 'source')
        outdir = os.path.join(os.path.abspath(self.build_lib),
                              'spacepy', 'Doc')
        cmd = '{0} -b html -d {1} {2} {3}'.format(
            os.environ['SPHINXBUILD'] if 'SPHINXBUILD' in os.environ
            else 'sphinx-build',
            builddir, indir, outdir)
        try:
            subprocess.check_call(cmd.split())
        except:
            warnings.warn(
                "Building docs failed. Help will not be available.")

    def run(self):
        """Actually perform the build"""
        self.compile_libspacepy()
        if sys.platform == 'win32':
            #Copy mingw32 DLLs. This keeps them around if ming is uninstalled,
            #but more important puts them where bdist_wininst and bdist_wheel
            #will include them in binary installers
            copy_dlls(os.path.join(self.build_lib, 'spacepy', 'mingw'))
        _build.run(self) #need subcommands BEFORE building irbem
        self.compile_irbempy()
        delete_old_files(self.build_lib)
        if self.build_docs:
            self.make_docs()
        else:
            self.copy_docs()


class install(_install):
    """Extends base distutils install to fix compiler options"""

    user_options = _install.user_options + compiler_options + [
        ('build-docs', None,
         'Build documentation with Sphinx (default: copy pre-built) [False]'),
        ]

    def initialize_options(self):
        self.build_docs = False
        initialize_compiler_options(self)
        _install.initialize_options(self)

    def finalize_options(self):
        #Because we are building extension modules ourselves, distutils
        #can't tell this is non-pure and we need to use platlib.
        if self.install_lib is None:
            self.install_lib = self.install_platlib
        _install.finalize_options(self)
        finalize_compiler_options(self)

    def get_outputs(self):
        """Tell distutils about files we put in build by hand"""
        outputs = _install.get_outputs(self)
        docs = [
            os.path.join(
                self.install_libbase, dirpath[len(self.build_lib) + 1:], f)
            for (dirpath, dirnames, filenames)
            in os.walk(os.path.join(self.build_lib, 'spacepy', 'Doc'))
            for f in filenames]
        #This is just so we know what a shared library is called
        comp = distutils.ccompiler.new_compiler(compiler=self.compiler)
        if hasattr(distutils.ccompiler, 'customize_compiler'):
            distutils.ccompiler.customize_compiler(comp)
        else:
            distutils.sysconfig.customize_compiler(comp)
        libspacepy = os.path.join(
            'spacepy', comp.library_filename('spacepy', lib_type='shared'))
        if os.path.exists(os.path.join(self.build_lib, libspacepy)):
            spacepylibs = [os.path.join(self.install_libbase, libspacepy)]
        else:
            spacepylibs = []
        irbemlibfiles = [os.path.join('spacepy', 'irbempy', f)
                         for f in get_irbem_libfiles()]
        irbemlibs = [
            os.path.join(self.install_libbase, f) for f in irbemlibfiles
            if os.path.exists(os.path.join(self.build_lib, f))]
        return outputs + docs + spacepylibs + irbemlibs


def copy_dlls(outdir):
    """Copy the mingw runtime libraries into a build

    :param str outdir: Final target directory of the DLLs in the build.
    """
    libdir = None
    libnames = None
    libneeded = ('libgfortran', 'libgcc_s', 'libquadmath', 'libwinpthread')
    for p in os.environ['PATH'].split(';'):
        if not os.path.isdir(p):
            continue
        libnames = [
            f for f in os.listdir(p) if f[-4:].lower() == '.dll'
            and f.startswith(libneeded)]
        if len(libnames) == len(libneeded):
            libdir = p
            break
    if libdir is None:
        raise RuntimeError("Can't locate runtime libraries.")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for f in libnames:
        shutil.copy(os.path.join(libdir, f), outdir)


class bdist_wininst(_bdist_wininst):
    """Handle compiler options, libraries for build on Windows install"""

    user_options = _bdist_wininst.user_options + compiler_options

    def initialize_options(self):
        initialize_compiler_options(self)
        _bdist_wininst.initialize_options(self)

    def finalize_options(self):
        _bdist_wininst.finalize_options(self)
        finalize_compiler_options(self)


if 'bdist_wheel' in sys.argv:
    class bdist_wheel(_bdist_wheel):
        """Handle compiler options for wheel build on Windows"""

        user_options = _bdist_wheel.user_options + compiler_options

        def initialize_options(self):
            initialize_compiler_options(self)
            _bdist_wheel.initialize_options(self)

        def finalize_options(self):
            _bdist_wheel.finalize_options(self)
            finalize_compiler_options(self)
            #Force platform-specific build (wheel finalize does this based
            #on explicitly declaring extension modules, and we handle them
            #by hand)
            #TODO: Complains config variable Py_DEBUG and WITH_PYMALLOC are
            #unset; python ABI tag may be incorrect
            #https://github.com/pypa/pip/issues/3383
            self.root_is_pure = False


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
        rebuild_static_docs()
        _sdist.run(self)


try:
    class config_fc(numpy.distutils.command.config_compiler.config_fc):
        """Get the options sharing with build, install"""

        def initialize_options(self):
            initialize_compiler_options(self)
            numpy.distutils.command.config_compiler.config_fc.initialize_options(
                self)

        def finalize_options(self):
            numpy.distutils.command.config_compiler.config_fc.finalize_options(
                self)
            finalize_compiler_options(self)
except:
    if egginfo_only:
        pass
    else:
        raise


packages = ['spacepy', 'spacepy.irbempy', 'spacepy.pycdf',
            'spacepy.plot', 'spacepy.pybats', 'spacepy.toolbox', ]
#If adding to package_data, also put in MANIFEST.in
package_data = ['data/*.*', 'pybats/sample_data/*', 'data/LANLstar/*', 'data/TS07D/TAIL_PAR/*']

setup_kwargs = {
    'name': 'spacepy',
    'version': '0.2.2',
    'description': 'SpacePy: Tools for Space Science Applications',
    'long_description': 'SpacePy: Tools for Space Science Applications',
    'author': 'SpacePy team',
    'author_email': 'spacepy@lanl.gov',
    'maintainer': 'Steve Morley, Josef Koller, Dan Welling, Brian Larsen, Mike Henderson, Jon Niehof',
    'maintainer_email': 'spacepy@lanl.gov',
    'url': 'https://github.com/spacepy/spacepy',
#download_url will override pypi, so leave it out http://stackoverflow.com/questions/17627343/why-is-my-package-not-pulling-download-url
#    'download_url': 'https://sourceforge.net/projects/spacepy/files/spacepy/',
    'requires': ['numpy (>=1.6)', 'scipy (>=0.10)', 'matplotlib (>=1.5)', 'python_dateutil',
                 'h5py', 'python (>=2.7, !=3.0)'],
    'packages': packages,
    'package_data': {'spacepy': package_data},
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
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    'keywords': ['magnetosphere', 'plasma', 'physics', 'space', 'solar.wind', 'space.weather', 'magnetohydrodynamics'],
    'license':  'PSF',
    'platforms':  ['Windows', 'Linux', 'MacOS X', 'Unix'],
    'cmdclass': {'build': build,
                 'install': install,
                 'bdist_wininst': bdist_wininst,
                 'sdist': sdist,
          },
}

if not egginfo_only:
    setup_kwargs['cmdclass']['config_fc'] = config_fc

if use_setuptools:
#Sadly the format here is DIFFERENT than the distutils format
    setup_kwargs['install_requires'] = [
        'numpy>=1.10,!=1.15.0',
        'scipy>=0.11',
        'matplotlib>=1.5',
        'h5py>=2.6',
        #Do not install ffnet on Windows since there's no binary
        #(people must hand-install)
        'ffnet>=0.7;platform_system!="Windows"',
        #ffnet needs networkx but not marked as requires, so to get it via pip
        #we need to ask for it ourselves
        'networkx>=1.0',
        'python_dateutil>=1.4',
        # AstroPy is only required to convert to/from AstroPy, so either
        # user has it or don't care.
        #'astropy>=1.0',
    ]
if 'bdist_wheel' in sys.argv:
    setup_kwargs['cmdclass']['bdist_wheel'] = bdist_wheel

# run setup from distutil
with warnings.catch_warnings(record=True) as warnlist:
    setup(**setup_kwargs)
    if len(warnlist):
        print('\nsetup produced the following warnings. '
              'Some functionality may be missing.\n')
        for w in warnlist:
            print(w.message)
