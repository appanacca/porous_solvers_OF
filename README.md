# Porous OpenFOAM miscellaneous

The repository contains some of the OpenFOAM extension that I developped for my PhD.
The naming convention is the same one used in the actual [manuscript](https://github.com/appanacca/manuscript).
It mainly contain solvers and utilities that are needed to study fibrous porous media. It use a different approach that the openFoam porous capability.

## Table of content

Each one of the different directory contains respectively:
+ **solvers**: contains the solver that implements the various equations for the Volume Average Method
+ **utilities**: contain two different pre-processing utilities that helps to initialize the porosity and the permeability tensor for each case
+ **cases**: the openaFoam completa cases that use the previous two developped executables
+ **database_kriging**: contain the python scripts and the input files that are needed to build the kriging metamodel for the effective permeability

## Installation
All the solvers run on OpenFOAM 3.0.
Some of the test cases initialization need also PyFOAM 0.6.6.

Each solver and utility needs to be installed following the classical openFoam pattern:

```bash
>> cd your_solver_directory/
>> wmake
```

## Usage
The esage of each different executable is documented in the wiki.