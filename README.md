# Porous OpenFOAM solver and utils

The repository contains some of the OpenFOAM extension that I developped for my PhD.
The naming convention is the same one used in the actual [manuscript](https://github.com/appanacca/manuscript).


## ico_FFOAM

I a solver for the effective permeability tensor in porous media, based on the icoFoam solver.
It represent a correction to the "classical" permeability tensor in presence of inertia.
The solver solves one row of the tensor, and the parameter **delta_ij** decide wich one.
For example if: delta_ij = (1 0 0) the solver output are the F_11, F_12, F13 components of the tensor.
This parameter needs to be specified in the __constant/transportProperties__ dictionary of the case.

## ico_KFOAM

Is the same as the solver above but now the output tensor is the permeability tensor in the laminar regime without inertia.

## VANS_variablePor

Is a solver based on icoFOAM that solves for a coupled free flow and porous media.
The porous field **por** and the effective permeability field **K** should be defined before run the solver (as usual).
The equation solved are the VANS equation for rigid porous media with varying porosity.
This solver actually couple the Navier-Stokes system and the VANS system using a penalization method.

Two other parameters should be defined in the __constant/transportProperties__ 
+ **Len_ref** is the reference lenght of the REV, it's used to compute the pore Reynolds number of each porous REV
+ **flag_k_update** is a flag to activate the metamodel for the permeability tensor. The metamodel search in the database file **result.dat** that should be present in the case folder 

The database search of the metamodel is implemented using the awsome [**nanoflann**](https://github.com/jlblancoc/nanoflann) library. It is a great library that implements efficent Nearest Neighbor search wih KD-trees data structures in C++. Feel free to explore it.

## set_por_hill_y and set_K_hill_y 

Are two utilities used to initialize the **K** and **por** fields by the user.
I know that the swak4foam library do the same job but for some reason I am not confortable at using it so I write these two simple utilities.
Also it could be modularized a little more, for example linking a function that define the interface instead that hardcoding it inside the utility. Feel free to improve my work.