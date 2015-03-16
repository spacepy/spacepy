#!/usr/bin/env python

"""Build and test SpacePy under Wine."""

import glob
import os
import os.path
import shutil
import subprocess
import sys
import tempfile

#https://github.com/certik/numpy-vendor/blob/master/setup-wine.sh
#http://technet.microsoft.com/en-us/library/cc759262%28v=ws.10%29.aspx

class Winerunner(object):
    def __init__(self):
        self.tmpdir = tempfile.mkdtemp()
        self.env = {}
        #Consider cutting this WAAAAY back for security reasons
        self.env[64] = os.environ.copy()
#Uncomment this once we have completely automated the numpy install
#        if 'DISPLAY' in self.env[64]:
#            del self.env[64]['DISPLAY']
        self.env[32] = self.env[64].copy()
        self.env[64]['WINEPREFIX'] = os.path.join(self.tmpdir, 'wine64')
        self.env[64]['WINEARCH'] = 'win64'
        self.env[32]['WINEPREFIX'] = os.path.join(self.tmpdir, 'wine32')
        self.env[32]['WINEARCH'] = 'win32'

        self.logfd = os.open(os.path.join(self.tmpdir, 'log.txt'),
                        os.O_WRONLY | os.O_CREAT)
        self.run(['wineboot'])

    def run(self, cmdargs, bits=[32, 64], **kwargs):
        try:
            iter(bits)
        except:
            bits = [bits]
            #check_call is pointless, almost always something is slightly wrong
        return [subprocess.call(['wine'] + cmdargs, env=self.env[b],
                                stdout=self.logfd, stderr=self.logfd, **kwargs)
                for b in bits]

    def finish(self, clean=True):
        os.close(self.logfd)
        if clean:
            shutil.rmtree(self.tmpdir)
        else:
            print(self.tmpdir)

    def c_drive(self, bits=32):
        return os.path.join(self.tmpdir, 'wine{0}'.format(bits), 'drive_c')

    def reginsert(self, contents, bits=32):
        p = subprocess.Popen(['wine', 'regedit', '-'], env=self.env[bits],
                             stdout=self.logfd, stderr=self.logfd,
                             stdin=subprocess.PIPE)
        p.communicate('REGEDIT4\n\n' + contents)
        p.wait()

installers = '/home/jnik/work/spacepy_win/installers'
if not os.path.exists(installers):
    installers = '/home/jtniehof/spacepy_win/installers'
runner = Winerunner()
#ALLUSERS doesn't appear to work on 64-bit
runner.run(['msiexec', '/i', os.path.join(installers, 'python-2.7.9.amd64.msi'),
           '/qn', 'TARGETDIR=C:\Python27'], 64)
runner.run(['msiexec', '/i', os.path.join(installers, 'python-2.7.9.msi'),
           '/qn', 'TARGETDIR=C:\Python27', 'ALLUSERS=1'], 32)

#ah, numpy installers...
#https://github.com/numpy/numpy/issues/572
#https://github.com/numpy/numpy/pull/4727
#https://github.com/scikit-learn/scikit-learn/pull/3256
#http://stackoverflow.com/questions/10167302/install-numpy-exe-silently-with-batch-file-in-windows
#numpy and scipy are superpacks and have two levels of "can't make quiet"
#others are simply bdist_wininst, which have no quiet support
runner.run([os.path.join(installers,
                         'numpy-1.7.2-win32-superpack-python2.7.exe'), '/s'],
           32)
runner.run([os.path.join(installers,
                         'scipy-0.12.1-win32-superpack-python2.7.exe'), '/s'],
           32)
#mpl 1.4.2, 1.3.1 doesn't include dependencies
runner.run([os.path.join(installers,
                         'matplotlib-1.2.1.win32-py2.7.exe'), '/s'],
           32)
runner.run([os.path.join(installers,
                         'h5py-2.4.0.win32-py2.7.exe'), '/s'],
           32)
#networkx 1.9.x requires decorator module
runner.run([os.path.join(installers,
                         'networkx-1.8.1-1_py27.exe'), '/s'],
           32) #not bdist_wininst, but silent doesn't work
runner.run([os.path.join(installers,
                         'ffnet-0.7.1.win32-py2.7.exe'), '/s'],
           32)
#THIS at least can go unattended
runner.run([os.path.join(installers,
                         'cdf35_0_2-setup-32.exe'), '/q'],
           32)

#oh, dear, mingw
#http://www.mingw.org/wiki/InstallationHOWTOforMinGW
#but don't use links in that guide!
#Despite what it says, maybe need dev for gmp, mpc, mpfr, zlib?
#binutils-2.24-1-mingw32-bin.tar.xz
#gcc-core-4.8.1-4-mingw32-bin.tar.lzma
#gcc-core-4.8.1-4-mingw32-dev.tar.lzma
#gcc-core-4.8.1-4-mingw32-dll.tar.lzma
#gcc-fortran-4.8.1-4-mingw32-bin.tar.lzma
#gcc-fortran-4.8.1-4-mingw32-dev.tar.lzma
#gcc-fortran-4.8.1-4-mingw32-dll.tar.lzma
#gettext-0.18.3.2-1-1-mingw32-bin.tar.xz
#gmp-5.1.2-1-mingw32-dev.tar.lzma
#gmp-5.1.2-1-mingw32-dll.tar.lzma
#libasprintf-0.18.3.2-1-mingw32-dll-0.tar.xz
#libgettextpo-0.18.3.2-1-mingw32-dll-0.tar.xz
#libiconv-1.14-3-mingw32-dll.tar.lzma
#libintl-0.18.3.2-1-mingw32-dll-8.tar.xz
#mingwrt-3.21-mingw32-dev.tar.xz
#mingwrt-3.21-mingw32-dll.tar.xz
#mpc-1.0.1-2-mingw32-dev.tar.lzma
#mpc-1.0.1-2-mingw32-dll.tar.lzma
#mpfr-3.1.2-2-mingw32-dev.tar.lzma
#mpfr-3.1.2-2-mingw32-dll.tar.lzma
#pthreads-w32-2.9.1-1-mingw32-dev.tar.lzma
#pthreads-w32-2.9.1-1-mingw32-dll.tar.lzma
#w32api-3.17-2-mingw32-dev.tar.lzma
#zlib-1.2.8-1-mingw32-dev.tar.lzma
#zlib-1.2.8-1-mingw32-dll.tar.lzma

c = runner.c_drive(32)
mingout = os.path.join(c, 'MinGW')
os.mkdir(mingout)
mingin = os.path.join(installers, 'mingw')
for f in os.listdir(mingin):
    cmd = ['tar', '-C', mingout, '-xf', os.path.join(mingin, f)]
    subprocess.check_call(cmd, stdout=runner.logfd, stderr=runner.logfd)
#Putting in python scripts dir even though direct call to f2py fails
#also apparently the case matters sometimes.
runner.reginsert('[HKEY_CURRENT_USER\Environment]\n'
                 '"PATH"="C:\\\\MinGW\\\\bin;C:\\\\Python27;'
                 'C:\\\\Python27\\\\Scripts"\n')
runner.reginsert('[HKEY_CURRENT_USER\Environment]\n'
                 '"Path"="C:\\\\MinGW\\\\bin;C:\\\\Python27;'
                 'C:\\\\Python27\\\\Scripts"\n')
#f2py dies horribly in link if dll's aren't in same place as cc1.exe, as.exe
#even though same gcc command line works fine if run by hand.
for d in ('libgcc_s_dw2-1.dll', 'libgmp-10.dll', 'libmpc-3.dll',
          'libmpfr-4.dll', 'libintl-8.dll', 'libiconv-2.dll',
          'zlib1.dll'):
    shutil.copy(os.path.join(mingout, 'bin', d),
                os.path.join(mingout, 'libexec', 'gcc', 'mingw32',
                             '4.8.1'))
    shutil.copy(os.path.join(mingout, 'bin', d),
                os.path.join(mingout, 'mingw32', 'bin'))

spacepybuild = os.path.join(c, 'spacepy')
spzip = sorted(glob.glob(os.path.join(installers, 'spacepy-*.zip')))[-1]
subprocess.check_call(['unzip', spzip, '-d', spacepybuild],
                      stdout=runner.logfd, stderr=runner.logfd)
spacepybuild = os.path.join(spacepybuild,
                            os.path.basename(spzip)[:-4])
#Need to specify a compiler, or do the distutils.cfg?
#The f2py thing is weird. The script being on the path is
#insufficient for the call in wine, but seemed to work on the
#XP VM. This should be fixed in the installer.
runner.run(['python', 'setup.py', 'bdist_wininst',
            '--f2py=python C:\\\\Python27\\\\Scripts\\\\f2py.py'],
           bits=32, cwd=spacepybuild)
#Now install it, grab the omni database, run unit tests
spinstaller = os.path.basename(
    sorted(glob.glob(os.path.join(spacepybuild, 'dist', 'spacepy-*.exe')))[-1])
runner.run([spinstaller, '/s'], bits=32,
           cwd=os.path.join(spacepybuild, 'dist'))
#Force the directory for the ".spacepy" file, so can populate it w/OMNI, etc.
runner.env[32]['SPACEPY'] = "C:\\spacepy_user"
dotspacepy = os.path.join(c, 'spacepy_user', '.spacepy')
os.makedirs(dotspacepy)
import spacepy #this is in the BUILD environment, not the target!
shutil.copytree(os.path.join(spacepy.DOT_FLN, 'data'),
                os.path.join(dotspacepy, 'data'))
retval = runner.run(['python', 'test_spacepy.py'], bits=32,
                    cwd=os.path.join(spacepybuild, 'tests'))[0]
if retval:
    print('Unit tests failed.')
    runner.finish(False)
else:
    print('{0} built and tested.'.format(spinstaller))
    os.copy(os.path.join(spacepybuild, 'dist', spinstaller),
            installers)
    runner.finish(True)
