#!/bin/sh
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=gsienkf
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=20
#SBATCH --job-name="stoch_unit_tests"

source ./module-setup.sh
module purge
module use $( pwd -P )
module load modules.stoch
EXEC=standalone_ca.x
# compile codes
sh compile_standalone_ca.hera_intel
if [ ! -f $EXEC ];then
  echo "compilation errors"
  exit 1
fi
#sh compile_compare_ca.sh

# copy input directory
#cp -r /scratch2/BMC/gsienkf/Philip.Pegion/stochastic_physics_unit_tests/input_data INPUT

#layout 1x1
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/1/g" input.nml
sed -i -e "s/NPY/1/g" input.nml
sed -i -e "s/CA_SGS/.true./g" input.nml
sed -i -e "s/CA_GLOBAL/.false/g" input.nml
sed -i -e "s/_STOCHINI_/.false./g" input.nml
export OMP_NUM_THREADS=1
time srun --label -n 6 $EXEC >& stdout.1x1
mkdir ca_layout_1x1
mv ca_out* ca_layout_1x1
# test 3 different domain decompositions and compare to baseline
#layout 1x4
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/1/g" input.nml
sed -i -e "s/NPY/4/g" input.nml
sed -i -e "s/CA_SGS/.true./g" input.nml
sed -i -e "s/CA_GLOBAL/.false/g" input.nml
sed -i -e "s/_STOCHINI_/.false./g" input.nml
export OMP_NUM_THREADS=1
time srun --label -n 24 $EXEC  >& stdout.1x4
mkdir ca_layout_1x4
mv ca_out* ca_layout_1x4
exit
#layout 2x2
export OMP_NUM_THREADS=2
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/2/g" input.nml
sed -i -e "s/NPY/2/g" input.nml
sed -i -e "s/CA_SGS/.true./g" input.nml
sed -i -e "s/CA_GLOBAL/.false/g" input.nml
sed -i -e "s/_STOCHINI_/.false./g" input.nml
time srun -n 24 $EXEC
mkdir ca_layout_2x2
mv ca_out* ca_layout_2x2
#layout 1x4
export OMP_NUM_THREADS=1
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/4/g" input.nml
sed -i -e "s/NPY/1/g" input.nml
sed -i -e "s/CA_SGS/.true./g" input.nml
sed -i -e "s/CA_GLOBAL/.false/g" input.nml
sed -i -e "s/_STOCHINI_/.false./g" input.nml
time srun -n 24 $EXEC
mkdir ca_layout_4x1
mv ca_out* ca_layout_4x1
# restart run
mv stochy_middle.nc INPUT/atm_stoch.res.nc
export OMP_NUM_THREADS=2
cp input.nml.ca_template input.nml
sed -i -e "s/NPX/4/g" input.nml
sed -i -e "s/NPY/1/g" input.nml
sed -i -e "s/CA_SGS/.true./g" input.nml
sed -i -e "s/CA_GLOBAL/.false/g" input.nml
sed -i -e "s/_STOCHINI_/.true./g" input.nml
time srun -n 24 $EXEC
rm ca_out* 

compare_output
if [ $? -ne 0 ];then
   echo "unit tests failed"
else
#   diff stochy_final.nc stochy_final_2.nc
   if [ $? -eq 0 ];then
      echo "unit tests successful"
      rm -rf ca_layout_*
      rm logfile*
      rm stochy*nc
      rm ../*.o ../*.mod
      rm ../libstochastic_physics.a
      rm $EXEC
   else
      echo "restart test failed"
   fi
fi
