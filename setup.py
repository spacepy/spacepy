#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# setup.py to install spacepy

__version__ = "$Revision: 1.21 $, $Date: 2010/08/05 22:38:51 $"
__author__ = 'The SpacePy Team, Los Alamos National Lab (spacepy@lanl.gov)'

# -------------------------------------
def compile_pybats():
     import os

     os.system('f2py -c spacepy/pybats/ctrace2d.pyf spacepy/pybats/trace2d.c')

def compile_oneralib():
	
	# compile oneralib
	import os, sys
	
	
	os.chdir('spacepy/onerapy/onera_lib_V4.1')
	
	F90files = ['source/onera_desp_lib.f', 'source/CoordTrans.f']
	functions = ['make_lstar1', 'make_lstar_shell_splitting1', \
		'coord_trans1 find_magequator1', 'find_mirror_point1', 
		'get_field1']
	
	# call f2py
	os.system('f2py --overwrite-signature -m onerapylib -h onerapylib.pyf '+' '.join(F90files) \
			+' only: ' + ' '.join(functions) + ' :')
	
	# intent(out) substitute list
	outlist = ['lm', 'lstar', 'blocal', 'bmin', 'xj', 'mlt', 'xout', 'bmin', 'posit', \
			'xgeo', 'bmir', 'bl', 'bxgeo']
	
	inlist = ['sysaxesin', 'sysaxesout', 'iyr', 'idoy', 'secs', 'xin']
	
	fln = 'onerapylib.pyf'
	
	print 'Substituting fortran intent(in/out) statements'
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
   		os.system('gfortran -c -w -O2 -fPIC *.f')
   		os.system('libtool -static -o libBL2.a *.o')
   		os.chdir('..')
   		os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')
   		
	elif sys.platform == 'linux2': # then linux   	
		os.system('gfortran -c -w -O2 -fPIC *.f')
		os.system('ar -r libBL2.a *.o')
		os.system('ranlib libBL2.a')
		os.chdir('..')
		os.system('f2py -c onerapylib.pyf source/onera_desp_lib.f -Lsource -lBL2 --fcompiler=gnu95')
	else:
		print sys.platform, ' not supported at this time'
		sys.exit(1)
	
	os.system('mv -f onerapylib.so ../')
	os.chdir('../../..')
	
	return

# -------------------------------------
def subst(pattern, replacement, filestr,
          pattern_matching_modifiers=None):
          
    """
    replace pattern by replacement in file
    pattern_matching_modifiers: re.DOTALL, re.MULTILINE, etc.
    """
    
    import os, re, sys, shutil
    
    
    if pattern_matching_modifiers is not None:
    	cp = re.compile(pattern, pattern_matching_modifiers)
    else:
    	cp = re.compile(pattern)
	
	if cp.search(filestr):  # any occurence of pattern?
		filestr = cp.sub(replacement, filestr)
		
    return filestr
	
# -------------------------------------
from distutils.core import setup
import os, sys, shutil

#import tooblox by reading file from repository
# this will provide mostly the query_yes_no function
execfile('spacepy/toolbox.py')

#test for python version 2.x where x>=5
try:
    dum = sys.version_info
    assert dum[0]==2
    assert dum[1]>=5
    import numpy
    import scipy
    import matplotlib
except:
    raise Exception("""SpacePy requires Python 2.X, where X>=5.\n
    Numpy, Scipy and Matplotlib(>=0.99) are also required\n
    Please install suitable versions.""")

# run compile for onera_desp_lib first
if os.path.exists('spacepy/onerapy/onerapylib.so'):
    ans = query_yes_no('\nDo you want to recompile the ONERA-DESP library?', default="no")
    if ans=='yes':
        compile_oneralib()
else:
    compile_oneralib()

# Compile PyBats
compile_pybats()

# create .spacepy in $HOME and move data
# read-in .rc file first
execfile('spacepy/data/spacepy.rc')
if DOT_FLN[:2] == '~/': DOT_FLN = os.environ['HOME']+'/'+DOT_FLN[2:]
if DOT_FLN[-1] == '/': DOT_FLN = DOT_FLN[:-1]


if os.path.exists(DOT_FLN):
    ans = query_yes_no('\n'+DOT_FLN+' already exists. Do you want to start fresh?', default="no")
    if ans=='no':
        fresh_install = False
    else:
        fresh_install = True
        i = 0
        while 1:
            if os.path.exists(DOT_FLN+'.bak.'+str(i)):
                i = i+1
            else:
                shutil.move(DOT_FLN, DOT_FLN+'.bak.'+str(i))
                break
else:
    fresh_install = True

if fresh_install:
    os.mkdir(DOT_FLN)
    os.chmod(DOT_FLN, 0777)
    os.mkdir(DOT_FLN+'/data')
    os.chmod(DOT_FLN+'/data', 0777)
    shutil.copy('spacepy/data/spacepy.rc', DOT_FLN+'/')
    shutil.copy('spacepy/data/tai-utc.dat', DOT_FLN+'/data')

pkg_files = ['onerapy/onerapylib.so','onerapy/*.py', 'doc/*.*', 'pybats/*.py', 'pybats/*.so']

# run setup from distutil
setup(name='spacepy',
      version='0.1',
      description='SpacePy: Tools for Space Science Applications',
      author='Steve Morley, Josef Koller, Dan Welling, Brian Larsen, Mike Henderson',
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
          ]
      )

# update/download packages
if fresh_install:
    dir = update()
    print("Data installed to " + dir)
else:
    ans = query_yes_no("\nDo you want to update OMNI database and leap seconds table? (Internet connection required)", default = "no")
    if ans=='yes':
        dir = update()
    else:
        print "\nRemember to update OMNI and leap seconds table occasionally by running spacepy.toolbox.update()"


# offer testing routine
ans = query_yes_no("\nDo you want to test your SpacePy installation", default="no")
if ans=='yes':
    import spacepy
    spacepy.test_all()

print "\nThanks for installing SpacePy."
