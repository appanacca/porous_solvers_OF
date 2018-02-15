cp -r  0_org/ 0
transformPoints -yawPitchRoll '(0 0 -90)'
blockMesh
set_por
set_K
tar -xvzf kriging_K.dat.tar.gz
VANS_variabPor
