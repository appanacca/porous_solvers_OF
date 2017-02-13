rm 0/p.template 0/U.template
decomposePar 
mpirun -n 20 renumberMesh -overwrite -parallel
