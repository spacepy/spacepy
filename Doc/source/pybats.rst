
########################################
PyBats - SWMF & BATS-R-US Analysis Tools
########################################

.. currentmodule:: spacepy.pybats

.. automodule:: spacepy.pybats

###################
Subpackage Contents
###################
		
.. rubric:: Submodules
There are submodules for most models included within the SWMF.  The classes
and methods contained within are code-specific, yielding power and
convenience at the cost of flexibility.

.. autosummary::
   :template: clean_module.rst
   :toctree: autosummary

   bats
   dgcpm
   dipole
..   gitm
   kyoto
   pwom
   ram
   rim
   
.. rubric:: Generic Classes
Top-level PyBats classes handle common-format input and output from the SWMF
and are very flexible.  However, they do little beyond open files for the user.
	    
.. autosummary::
   :template: clean_class.rst
   :toctree: autosummary

   IdlBin
   ImfInput
   LogFile
   NgdcIndex
   PbData
   SatOrbit
