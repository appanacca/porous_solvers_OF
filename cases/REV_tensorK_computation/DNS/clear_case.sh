#!/bin/bash

pyFoamClearCase.py . --function-object-data --clear-history --remove-analyzed --processors-remove
foamClearPolyMesh
rm Ply*
rm PyFoamHistory
rm *.log
rm *.png
rm *.dat
rm -r 0
shopt -s extglob
cd constant/polyMesh/
rm !(blockMeshDict.template)
cd ../..
rm caseSetup.sh
rm PyFoam*
rm meshCreate.sh
rm *.eps
