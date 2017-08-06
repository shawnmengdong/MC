# MC
MatrixCondensation

This is a matrix condensation repository developed by XM.Dong, all algorithms are in Fortran95

the source code includes:

mc_s : serial version of matirix condensation  (which is essentially same as gaussian elimination)
mc_p : parallel version of matrix condensation using MPI

ge_s : serial version of gaussian elimination
ge_p : parallel version of gaussian elimination using local partial pivoting 

ge_pg: parallel version of gaussian elimination using global partial pivoting


the ouput folder is the experiment run from OSCER cluster(OU Supercomputing Center for Education&Research)


How to use:

Make sure you have openMPI and gfortran installed before compilation.


Example:

make mc_p

mpirun -n 4 mc_p 1000     (this will calculate log(abs(det)) of random 1000x1000 dense matrix)

mpirun -n 4 mc_p 1000 /root/to/specific/test/matrix  (this will calculate a log(abs(det)) of specified matrix data)

test matrix can be downloaded here: http://morpheus.mcs.utulsa.edu/~papama/hpc/



if you have any question, send email to Xiaomeng.dong-1@ou.edu

if you want to contribute, let me know. 

