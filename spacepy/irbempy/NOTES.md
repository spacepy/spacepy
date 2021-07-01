# Notes on how to update the revision of IRBEMlib
First you'll need to check out the IRBEM git repository from
https://github.com/PRBEM/IRBEM
Some required files get built by make, so in the IRBEM repo:
```
make OS=linux64 ENV=gnu64 all
```

The `compare_irbem.py` script in `spacepy/developer/scripts` will list the
files that differ between the IRBEM repo and the current spacepy version.
It will also give line numbers with unicode characters that need replacing.

### If necessary, update myOwnMagField.f
SpacePy uses the `myOwnMagneticField.f` to implement the "Dungey" model.
In this case that's a centered dipole plus uniform field.
In the unlikely event that there's a required change to the calling syntax
or similar, this file will need updating. Otherwise, it's expected to be
different since the file in the IRBEM directory is a stub.

### fix comments; unicode
It seems that some comment lines confuse f2py.
E.g., line 690 (in SVN rev546) "      REAL*8 Bmirror  ! particle's mirror field strength"
The apostrophe breaks the compile, so anything like this needs to be removed.
Some files have unicode degree symbols that should be removed. `compare_irbem.py` should flag these.

### TS07D updates
ts07d.inc is required to build for SpacePy, but is not in the repo on update.
If this is not present, then run make as detailed above and add to SpacePy's
IRBEM source folder.

## Changing SpacePy's IRBEM source

For each file highlighted by `compare_irbem.py`, *except*
`myOwnMagneticField.f` replace the file in
`spacepy/spacepy/irbempy/irbem-lib-[whatever the latest rev is]`

If a line is flagged as having unicode characters, head to those lines and
replace. E.g. if it's supposed to be a degree symbol, replace with "deg"

### setup.py
Finally, in `setup.py` replace the path to SpacePy's IRBEM source code
following the suggestion from `compare_irbem.py`. You'll need to do this
using `git mv`
