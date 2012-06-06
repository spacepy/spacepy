"""Gets a list of classes defined in every python file and spits out
__all__ line for each"""

import glob
import subprocess

classdefs = {}
for fname in glob.glob('*.py'):
    cmd = ['grep', '-e', r'^class ', fname]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    (outdata, errdata) = p.communicate()
    classlist = ["'" + l.split(' ')[1].split('(')[0] + "'"
                 for l in outdata.split('\n')
                 if l]
    if classlist:
        print('__all__ = [{1}] ({0})'.format(fname, ', '.join(classlist)))
        for c in classlist:
            if not c in classdefs:
                classdefs[c] = [fname]
            else:
                classdefs[c].append(fname)

for k in classdefs:
    if len(classdefs[k]) > 1:
        print('{0} multiply defined ({1})'.format(k, ', '.join(classdefs[k])))

