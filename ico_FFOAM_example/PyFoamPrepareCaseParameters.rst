.. title:: t0_f0_g1006_n50
.. sectnum::
.. header:: t0_f0_g1006_n50
.. header:: Wed Nov  2 17:04:07 2016
.. footer:: ###Page### / ###Total###

Parameters set in case directory ``/tmpdir/nluminar/parametric_K_database/t0_f0_g1006_n50`` at *Wed Nov  2 17:04:07 2016*

.. contents::


***************
Parameter files
***************
Parameters read from files

1. ``parameters.basic``



************************
Parameters with defaults
************************

Automatically defined values
============================

Short name: ``System``

================== ======================================================== ========================================================
..                 default                                                  Value                                                   
================== ======================================================== ========================================================
casePath           "/tmpdir/nluminar/parametric_K_database/t0_f0_g1006_n50" "/tmpdir/nluminar/parametric_K_database/t0_f0_g1006_n50"
foamFork           openfoamplus                                             openfoamplus                                            
caseName           "t0_f0_g1006_n50"                                        "t0_f0_g1006_n50"                                       
foamVersion        3.0                                                      3.0                                                     
numberOfProcessors 1                                                        1                                                       
================== ======================================================== ========================================================

**********************
Unspecified parameters
**********************
If these parameters are actually used then specify them in ``default.parameters``


=========== ============
..          Value       
=========== ============
phi         0           
g           1006        
n_elem_mesh 50          
porosity    0.6073009183
nproc       20          
theta       0           
D2          off         
D3          on          
=========== ============

*****************
Derived Variables
*****************
Script with derived Parameters``./derivedParameters.py``

These values were added:

== =====================================================================================================
.. Value                                                                                                
== =====================================================================================================
fx 1006.0                                                                                               
fy 0.0                                                                                                  
fz 0.0                                                                                                  
np <module 'numpy' from '/users/p1540/nluminar/anaconda/lib/python2.7/site-packages/numpy/__init__.pyc'>
R  0.500000000001                                                                                       
== =====================================================================================================
The code of the script:

.. code:: python
  import numpy as np
  
  fx = g*np.cos(theta*np.pi/180)*np.cos(phi*np.pi/180)
  fy = g*np.sin(theta*np.pi/180)*np.cos(phi*np.pi/180)
  fz = g*np.sin(phi*np.pi/180)
  
  R = np.sqrt((4*(1-porosity))/(2*np.pi))
  

