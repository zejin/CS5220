#!/bin/sh

for nproc in {1,2,4,8}
do
file="strong_${nproc}.pbs"
echo "#!/bin/sh -l" >> $file
echo "#PBS -l nodes=1:ppn=24" >> $file
echo "#PBS -l walltime=0:30:00" >> $file
echo "#PBS -N LTM_mpi" >> $file
echo "#PBS -j oe" >> $file
echo "module load cs5220" >> $file
echo "cd \$PBS_O_WORKDIR" >> $file

echo "echo 'strong scaling'">> $file
echo "echo 'number of processors: $nproc'" >> $file

echo "mpirun -np $nproc LTM_mpi" >> $file

qsub $file
rm $file
done