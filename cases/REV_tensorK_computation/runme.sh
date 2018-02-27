cd DNS

	pyFoamPrepareCase.py . --no-mesh-create  --parameter-file=parameters.basic
	mpirun -n 20 icoGFoam -parallel
	reconstructPar -latestTime
	scale_U
	python getResults.py postProcessing/forces/0/force.dat postProcessing/avgU/0/cellSource.dat 0.03 > ../resuts.dat
	rm -r processor*

cd ../1_row

	cp ../DNS/0.05/scaleU 0/
	mv 0/scaleU 0/U
	mpirun -20 ico_FFoam -parallel
	reconstructPar -latestTime
	python getResults.py postProcessing/avgU/0/cellSource.dat  0.03 >> ../resuts.dat
	rm -r processor*

cd ../2_row

	cp ../DNS/0.05/scaleU 0/
	mv 0/scaleU 0/U
	mpirun -n 20 ico_FFoam -parallel
	reconstructPar -latestTime
	python getResults.py postProcessing/avgU/0/cellSource.dat 0.03 >> ../resuts.dat
	rm -r processor*

cd ../3_row

	cp ../DNS/0.05/scaleU 0/
	mv 0/scaleU 0/U
	mpirun -n 20 ico_FFoam -parallel
	reconstructPar -latestTime
	python getResults.py postProcessing/avgU/0/cellSource.dat 0.03 >> ../resuts.dat
	rm -r processor*
