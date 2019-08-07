# Notes on how to update the revision of IRBEMlib

### convert .../source/drift_bounce_orbit.f to ASCII
iconv -f ISO-8859-1 -t ASCII -c drift_bounce_orbit.f > drift_bounce_orbit_new.f
mv drift_bounce_orbit_new.f drift_bounce_orbit.f

### replace myOwnMagField.f with Dungey model
cp [previousRev]/source/myOwnMagField.f [newRev]/source/myOwnMagField.f

### fix comments; unicode
It seems that some comment lines confuse f2py.
E.g., line 690 (in rev546) "      REAL*8 Bmirror  ! particle's mirror field strength"
The apostrophe breaks the compile, so it needs to be removed.
Lines 331 and 3083 (in rev616) have unicode degree symbols that should be removed.
There are a lot of wanrings about lines with an exclamation as comment marker - for now we ignore these.

### TS07D updates
ts07d.inc is required to build for SpacePy, but is not in the repo on update.
Run make on SVN repo to generate the .inc file.

In SVN tree:
make OS=linux64 ENV=gnu64 all #sub as appropriate

then:
cp [IRBEM SVN tree]/source/ts07d.inc [newRev]/source

### setup.py
replace path to IRBEM source code
