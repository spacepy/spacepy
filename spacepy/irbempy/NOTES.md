# Notes on how to update the revision of IRBEMlib

### convert .../source/drift_bounce_orbit.f to ASCII
iconv -f ISO-8859-1 -t ASCII -c drift_bounce_orbit.f > drift_bounce_orbit_new.f
mv drift_bounce_orbit_new.f drift_bounce_orbit.f

NB: Not sure whether this is strictly necessary

### replace myOwnMagField.f with Dungey model
cp [previousRev]/source/myOwnMagField.f [newRev]/source/myOwnMagField.f

### fix comments
It seems that some comment lines confuse f2py.
E.g., line 689 (in rev541) "      REAL*8 Bmirror  ! particle's mirror field strength"
The apostrophe breaks the compile.
For safety I removed the apostrophe and converted blank lines starting with ! to start with c.

### TS07D updates
ts07d.inc is required to build for SpacePy, but is not in the repo on update.
Unpack the .tar.bz2 file to get the data, and run make on repo to generate the .inc.

In SVN tree:
bzip2 -dk TS07d.tar.bz2
tar -xvf TS07d.tar
make OS=linux64 ENV=gnu64 all #sub as appropriate

then:
cp [IRBEM SVN tree]/source/ts07d.inc [newRev]/source

### setup.py
replace path to IRBEM source code
